#!/usr/bin/env python
# coding: utf-8

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import seaborn as sns

class Reports:
    metadata_table = None
    frag_hist = None
    frag_violin = None
    frag_bin500 = None
    seacr_beds = None
    bams = None

    def __init__(self, logger, meta, meta_ctrl, raw_frags, bin_frag, seacr_bed):
        self.logger = logger
        self.meta_path = meta
        self.meta_ctrl_path = meta_ctrl
        self.raw_frag_path = raw_frags
        self.bin_frag_path = bin_frag
        self.seacr_bed_path = seacr_bed

        sns.set()
        sns.set_theme()
        sns.set_context("paper")

    #*
    #========================================================================================
    # UTIL
    #========================================================================================
    #*/

    def format_millions(self, x, pos):
        #the two args are the value and tick position
        return '%1.1fM' % (x * 1e-6)

    def format_thousands(self, x, pos):
        #the two args are the value and tick position
        return '%1.1fK' % (x * 1e-3)

    def gen_pdf(self, output_path, plots):
        with PdfPages(os.path.join(output_path, 'merged_report.pdf')) as pdf:
            for key in plots:
                pdf.savefig(plots[key])

    #*
    #========================================================================================
    # ENTRY POINTS
    #========================================================================================
    #*/

    def generate_cutandrun_reports(self, output_path):
        # Init
        abs_path = os.path.abspath(output_path)

        # Get plots and supporting data tables
        plots, data, mqc_frag_hist = self.generate_reports()

        # Save mqc text file
        txt_mqc = open(os.path.join(abs_path, "03_03_frag_len_mqc.txt"), "w")
        txt_mqc.write(mqc_frag_hist)
        txt_mqc.close()

        # Save data to output folder
        for key in data:
            data[key].to_csv(os.path.join(abs_path, key + '.csv'), index=False)

        # Save plots to output folder
        for key in plots:
            plots[key].savefig(os.path.join(abs_path, key + '.png'))

        # Save pdf of the plots
        self.gen_pdf(abs_path, plots)

    #*
    #========================================================================================
    # LOAD DATA
    #========================================================================================
    #*/

    def load_meta_data(self):
        # Plots supported

        # Sequencing Depth
        # Alignable Fragments
        # Alignment Rate (Target)
        # Alignment Rate (Spike-in)
        # Duplication Rate
        # Estimated Library Size
        # Unique Fragments
        # Spike-in Scale Factor
        # Normalised Fragment Count

        self.metadata_noctrl_table = pd.read_csv(self.meta_path, sep=',')
        self.metadata_table = pd.read_csv(self.meta_ctrl_path, sep=',')

        self.duplicate_info = False
        if 'dedup_percent_duplication' in self.metadata_table.columns:
            self.duplicate_info = True

        self.scale_factor_info = False
        if 'scale_factor' in self.metadata_table.columns:
            self.scale_factor_info = True

        # Make new perctenage alignment columns
        self.metadata_table['target_alignment_rate'] = self.metadata_table.loc[:, ('bt2_total_aligned_target')] / self.metadata_table.loc[:, ('bt2_total_reads_target')] * 100
        self.metadata_table['spikein_alignment_rate'] = self.metadata_table.loc[:, ('bt2_total_aligned_spikein')] / self.metadata_table.loc[:, ('bt2_total_reads_spikein')] * 100

        # Change to percentage
        if 'frip' in self.metadata_noctrl_table.columns:
            self.metadata_noctrl_table['frip'] =  self.metadata_noctrl_table['frip'] * 100

        # Sort tables
        self.metadata_table = self.metadata_table.sort_values('group')
        self.metadata_noctrl_table = self.metadata_noctrl_table.sort_values('group')

    def load_raw_frag_histogram(self):
        # Plots supported

        # Fragment Length Distribution (the main frag length lineplot)
        # Fragment Length Distribution (the main frag length violin plot)

        # Create list of deeptools raw fragment files
        dt_frag_list = glob.glob(self.raw_frag_path)
        dt_frag_list.sort()

        for i in list(range(len(dt_frag_list))):
            # Create dataframe from csv file for each file and save to a list
            dt_frag_i = pd.read_csv(dt_frag_list[i], sep='\t', header=None, names=['Size','Occurrences'])
            frag_base_i = os.path.basename(dt_frag_list[i])

            #  split txt files on dots
            sample_id_list = frag_base_i.split(".")

            # join list on the elements of the sample id
            separator = ""
            sample_id = separator.join(sample_id_list[0:-2])

            # split sample id on underscores
            sample_id_split_list = sample_id.split("_")

            #  take first element of this list for group id
            group_i = separator.join(sample_id_split_list[0:-1])

            #  take last element of this list for replicate number
            rep_i = sample_id_split_list[-1]

            # Create long forms of fragment histograms
            dt_frag_i_long = np.repeat(dt_frag_i['Size'].values, dt_frag_i['Occurrences'].values)
            dt_group_i_long = np.repeat(group_i, len(dt_frag_i_long))
            dt_rep_i_long = np.repeat(rep_i, len(dt_frag_i_long))

            dt_group_i_short = np.repeat(group_i, dt_frag_i.shape[0])
            dt_rep_i_short = np.repeat(rep_i, dt_frag_i.shape[0])

            if i==0:
                frags_arr = dt_frag_i_long
                group_arr = dt_group_i_long
                rep_arr = dt_rep_i_long

                group_short = dt_group_i_short
                rep_short = dt_rep_i_short
                self.frag_hist = dt_frag_i
            else:
                frags_arr = np.append(frags_arr, dt_frag_i_long)
                group_arr = np.append(group_arr, dt_group_i_long)
                rep_arr = np.append(rep_arr, dt_rep_i_long)

                group_short = np.append(group_short, dt_group_i_short)
                rep_short = np.append(rep_short, dt_rep_i_short)
                self.frag_hist = self.frag_hist.append(dt_frag_i)

        self.frag_hist['group'] = group_short
        self.frag_hist['replicate'] = rep_short
        self.frag_hist = self.frag_hist.reset_index(drop=True)
        self.frag_violin = pd.DataFrame( { "fragment_size" : frags_arr, "group" : group_arr , "replicate": rep_arr} )

    def load_binned_frags(self):
        # Plots supported

        # Replicate Reproducibility

        # Start by creating list of bin500 files
        dt_bin_frag_list = glob.glob(self.bin_frag_path)
        dt_bin_frag_list.sort()

        for i in list(range(len(dt_bin_frag_list))):
            dt_bin_frag_i_read = pd.read_csv(dt_bin_frag_list[i], sep='\t', header=None, names=['chrom','bin','count','sample'])
            sample_name = dt_bin_frag_i_read['sample'].iloc[0].split(".")[0]
            dt_bin_frag_i = dt_bin_frag_i_read[['chrom','bin','count']]
            dt_bin_frag_i.columns = ['chrom','bin',sample_name]

            if i==0:
                self.frag_bin500 = dt_bin_frag_i

            else:
                self.frag_bin500 = pd.merge(self.frag_bin500, dt_bin_frag_i, on=['chrom','bin'], how='outer')

        # Add log2 transformed count data column
        # log2_counts = self.frag_bin500[self.frag_bin500.columns[-(len(dt_bin_frag_list)):]].transform(lambda x: np.log2(x))
        # self.frag_bin500 = pd.concat([self.frag_bin500[['chrom','bin']],log2_counts], axis=1)

        self.frag_bin500 = self.frag_bin500.fillna(0)

    def load_seacr_peaks(self):
        # Plots supported

        # Total Peaks
        # Peak Width

        # Create dataframe for seacr peaks
        seacr_bed_list = glob.glob(self.seacr_bed_path)

        # combine all seacr bed files into one df including group and replicate info
        for i in list(range(len(seacr_bed_list))):
            seacr_bed_i = pd.read_csv(seacr_bed_list[i], sep='\t', header=None, usecols=[0,1,2,3,4], names=['chrom','start','end','total_signal','max_signal'])
            bed_base_i = os.path.basename(seacr_bed_list[i])

            #  split bed files on dots
            bed_id_list = bed_base_i.split(".")

            # join list on the elements of the sample id
            separator = ""
            sample_id = separator.join(bed_id_list[0:-4])

            # split sample id on underscores
            sample_id_split_list = sample_id.split("_")

            #  take first element of this list for group id
            group_i = separator.join(sample_id_split_list[0:-1])

            # take last element fo this list for replicate number
            rep_i = sample_id_split_list[-1]

            seacr_bed_i['group'] = np.repeat(group_i, seacr_bed_i.shape[0])
            seacr_bed_i['replicate'] = np.repeat(rep_i, seacr_bed_i.shape[0])

            if i==0:
                self.seacr_beds = seacr_bed_i

            else:
                self.seacr_beds = self.seacr_beds.append(seacr_bed_i)

        # Get peak stats by group and replicate
        self.seacr_beds_group_rep = self.seacr_beds[['group','replicate']].groupby(['group','replicate']).size().reset_index().rename(columns={0:'all_peaks'})

    def load_data(self):
        # ---------- Data - Load main meta-data table --------- #
        self.load_meta_data()

        # ---------- Data - Raw frag histogram --------- #
        self.load_raw_frag_histogram()

        # ---------- Data - Binned frags --------- #
        self.load_binned_frags()

        # ---------- Data - Peaks --------- #
        self.load_seacr_peaks()

    #*
    #========================================================================================
    # GEN REPORTS
    #========================================================================================
    #*/

    def generate_reports(self):
        # Init
        plots = dict()
        data = dict()

        # Get Data
        self.load_data()

        # Plot section 1
        multi_plot, data1 = self.alignment_summary()
        plots["01_01_seq_depth"] = multi_plot[0]
        plots["01_02_alignable_frag"] = multi_plot[1]
        plots["01_03_alignment_rate_target"] = multi_plot[2]
        plots["01_04_alignment_rate_spikein"] = multi_plot[3]
        data["01_alignment_summary"] = data1

        # Plot section 2
        if self.duplicate_info == True:
            multi_plot, data2 = self.duplication_summary()
            plots["02_01_dup_rate"] = multi_plot[0]
            plots["02_02_est_lib_size"] = multi_plot[1]
            plots["02_03_unique_frags"] = multi_plot[2]
            data["02_duplication_summary"] = data2

        #  Plot 3
        plot3, data3 = self.fraglen_summary_violin()
        plots["03_01_frag_len_violin"] = plot3
        data["03_01_frag_len_violin"] = data3

        # Plot 4
        plot4, data4 = self.fraglen_summary_histogram()
        plots["03_02_frag_len_hist"] = plot4
        data["03_02_frag_len_hist"] = data4

        # Plot 5
        plot5, data5 = self.replicate_heatmap()
        plots["04_replicate_heatmap"] = plot5
        data["04_replicate_heatmap"] = data5

        # Plot section 6
        if self.scale_factor_info == True:
            multi_plot, data6 = self.scale_factor_summary()
            plots["05_01_scale_factor"] = multi_plot[0]
            plots["05_02_frag_count"] = multi_plot[1]
            data["05_scale_factor_summary"] = data6

        # Plot 7a
        plot7a, data7a = self.no_of_peaks()
        plots["06_01_no_of_peaks"] = plot7a
        data["06_01_no_of_peaks"] = data7a

        # Plot 7b
        plot7b, data7b = self.peak_widths()
        plots["06_02_peak_widths"] = plot7b
        data["06_02_peak_widths"] = data7b

        # Plot 7c
        if 'peak_repro' in self.metadata_noctrl_table.columns:
            plot7c, data7c = self.reproduced_peaks()
            plots["06_03_reproduced_peaks"] = plot7c
            data["06_03_reproduced_peaks"] = data7c

        # Plot 7d
        if 'frip' in self.metadata_noctrl_table.columns:
            plot7d, data7d = self.frags_in_peaks()
            plots["06_04_frags_in_peaks"] = plot7d
            data["06_04_frags_in_peaks"] = data7d

        # Fragment Length Histogram data in MultiQC yaml format
        mqc_frag_hist = self.frag_len_hist_mqc()

        return (plots, data, mqc_frag_hist)

    #*
    #========================================================================================
    # TXT FILES
    #========================================================================================
    #*/

    def frag_len_hist_mqc(self):
        size_list = self.frag_hist["Size"].to_numpy().astype(str)
        occurrences_list = self.frag_hist["Occurrences"].to_numpy().astype(str)
        size_list_sep = np.core.defchararray.add(size_list, " : ")
        x_y_list = np.core.defchararray.add(size_list_sep, occurrences_list)

        group_rep = self.frag_hist[['group','replicate']].groupby(['group','replicate']).size().reset_index()
        first_line = "data:"

        for i in list(range(group_rep.shape[0])):
            group_i = group_rep.at[i,"group"]
            rep_i = group_rep.at[i,"replicate"]
            str_list = x_y_list[ (self.frag_hist['group'] == group_i) & (self.frag_hist['replicate'] == rep_i) ]
            # print(str_list)

            x_y_str = ", ".join(str_list)
            full_line_i = "    '" + group_i + "_" + rep_i + "' : {" + x_y_str + "}"
            if i==0:
                frag_len_hist_mqc_dict = "\n".join([first_line, full_line_i])

            else:
                frag_len_hist_mqc_dict = "\n".join([frag_len_hist_mqc_dict, full_line_i])


        return frag_len_hist_mqc_dict

    #*
    #========================================================================================
    # PLOTS
    #========================================================================================
    #*/

    # ---------- Plot 1 - Alignment Summary --------- #
    def alignment_summary(self):
        sns.color_palette("magma", as_cmap=True)
        sns.set(font_scale=0.6)
        # Subset data
        df_data = self.metadata_table.loc[:, ('id', 'group', 'bt2_total_reads_target', 'bt2_total_aligned_target', 'target_alignment_rate', 'spikein_alignment_rate')]

        # Create plots array
        figs = []

        # Seq depth
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_data, x='group', y='bt2_total_reads_target', palette = "magma")
        fig.suptitle("Sequencing Depth")
        ax.set(ylabel="Total Reads")
        ax.xaxis.set_tick_params(labelrotation=45)
        figs.append(fig)

        # Alignable fragments
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_data, x='group', y='bt2_total_aligned_target', palette = "magma")
        fig.suptitle("Alignable Fragments")
        ax.set(ylabel="Total Aligned Reads")
        ax.xaxis.set_tick_params(labelrotation=45)
        figs.append(fig)

        # Alignment rate hg38
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_data, x='group', y='target_alignment_rate', palette = "magma")
        fig.suptitle("Alignment Rate (Target)")
        ax.set(ylabel="Percent of Fragments Aligned")
        #ax.set(ylim=(0, 100))
        ax.xaxis.set_tick_params(labelrotation=45)
        figs.append(fig)

        # Alignment rate e.coli
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_data, x='group', y='spikein_alignment_rate', palette = "magma")
        fig.suptitle("Alignment Rate (Spike-in)")
        ax.set(ylabel="Percent of Fragments Aligned")
        #ax.set(ylim=(0, 100))
        ax.xaxis.set_tick_params(labelrotation=45)
        figs.append(fig)

        return figs, df_data

    # ---------- Plot 2 - Duplication Summary --------- #
    def duplication_summary(self):
        # Init
        k_formatter = FuncFormatter(self.format_thousands)
        m_formatter = FuncFormatter(self.format_millions)

        # Subset data
        df_data = self.metadata_table.loc[:, ('id', 'group', 'dedup_percent_duplication', 'dedup_estimated_library_size', 'dedup_read_pairs_examined')]
        df_data['dedup_percent_duplication'] *= 100
        df_data['unique_frag_num'] = df_data['dedup_read_pairs_examined'] * (1-df_data['dedup_percent_duplication']/100)

        # Create plots array
        figs = []

        # Duplication rate
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_data, x='group', y='dedup_percent_duplication', palette = "magma")
        fig.suptitle("Duplication Rate")
        ax.set(ylabel="Rate (%)")
        ax.set(ylim=(0, 100))
        ax.xaxis.set_tick_params(labelrotation=45)
        figs.append(fig)

        # Estimated library size
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_data, x='group', y='dedup_estimated_library_size', palette = "magma")
        fig.suptitle("Estimated Library Size")
        ax.set(ylabel="Size")
        ax.yaxis.set_major_formatter(m_formatter)
        ax.xaxis.set_tick_params(labelrotation=45)
        figs.append(fig)

        # No. of unique fragments
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_data, x='group', y='unique_frag_num', palette = "magma")
        fig.suptitle("Unique Fragments")
        ax.set(ylabel="Count")
        ax.yaxis.set_major_formatter(k_formatter)
        ax.xaxis.set_tick_params(labelrotation=45)
        ax.set(ylim=(0))
        figs.append(fig)

        return figs, df_data

    # ---------- Plot 3 - Fragment Distribution Violin --------- #
    def fraglen_summary_violin(self):
        fig, ax = plt.subplots()
        ax = sns.violinplot(data=self.frag_violin, x="group", y="fragment_size", hue="replicate", palette = "viridis")
        ax.set(ylabel="Fragment Size")
        ax.xaxis.set_tick_params(labelrotation=45)
        fig.suptitle("Fragment Length Distribution")

        return fig, self.frag_violin

    # ---------- Plot 4 - Fragment Distribution Histogram --------- #
    def fraglen_summary_histogram(self):
        fig, ax = plt.subplots()
        ax = sns.lineplot(data=self.frag_hist, x="Size", y="Occurrences", hue="group", style="replicate", palette = "magma")
        fig.suptitle("Fragment Length Distribution")

        return fig, self.frag_hist

    # ---------- Plot 5 - Replicate Reproducibility Heatmap --------- #
    def replicate_heatmap(self):
        fig, ax = plt.subplots()
        plot_data = self.frag_bin500[self.frag_bin500.columns[-(len(self.frag_bin500.columns)-2):]]
        # plot_data = plot_data.fillna(0)
        corr_mat = plot_data.corr(method='pearson')
        ax = sns.heatmap(corr_mat, annot=True, vmin=0, vmax=1)
        fig.suptitle("Replicate Reproducibility (read counts in 500bp bins)")

        return fig, self.frag_bin500

    # ---------- Plot 6 - Scale Factor Comparison --------- #
    def scale_factor_summary(self):
        # Get normalised count data
        df_normalised_frags = self.metadata_table.loc[:, ('id', 'group')]
        df_normalised_frags['normalised_frags'] = self.metadata_table['bt2_total_reads_target']*self.metadata_table['scale_factor']

        figs = []

        # Subset meta data
        df_data_scale = self.metadata_table.loc[:, ('id', 'group','scale_factor')]

        # Scale factor
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_data_scale, x='group', y='scale_factor', palette = "magma")
        fig.suptitle("Spike-in Scale Factor")
        ax.set(ylabel="Coefficient")
        ax.xaxis.set_tick_params(labelrotation=45)
        figs.append(fig)

        # Normalised fragment count
        fig, ax = plt.subplots()
        ax = sns.boxplot(data=df_normalised_frags, x='group', y='normalised_frags', palette = "magma")
        fig.suptitle("Normalised Fragment Count")
        ax.set(ylabel="Count")
        ax.xaxis.set_tick_params(labelrotation=45)
        figs.append(fig)

        return figs, df_data_scale

    # ---------- Plot 7 - Peak Analysis --------- #
    def no_of_peaks(self):
    # 7a - Number of peaks
        fig, ax = plt.subplots()
        fig.suptitle("Total Peaks")

        ax = sns.boxplot(data=self.seacr_beds_group_rep, x='group', y='all_peaks', palette = "magma")
        ax.xaxis.set_tick_params(labelrotation=45)
        ax.set_ylabel("No. of Peaks")

        return fig, self.seacr_beds_group_rep

    # 7b - Width of peaks
    def peak_widths(self):
        fig, ax = plt.subplots()

        ## add peak width column
        self.seacr_beds['peak_width'] = self.seacr_beds['end'] - self.seacr_beds['start']
        self.seacr_beds['peak_width'] = self.seacr_beds['peak_width'].abs().astype('float64')

        ax = sns.violinplot(data=self.seacr_beds, x="group", y="peak_width", hue="replicate", palette = "viridis")
        ax.xaxis.set_tick_params(labelrotation=45)
        ax.set_ylabel("Peak Width")
        fig.suptitle("Peak Width Distribution")

        return fig, self.seacr_beds


    # 7c - Peaks reproduced
    def reproduced_peaks(self):
        fig, ax = plt.subplots()

        # Subset data
        df_data = self.metadata_noctrl_table.loc[:, ('id', 'group', 'peak_repro')]

        # plot
        ax = sns.boxplot(data=df_data, x="group", y="peak_repro", palette = "magma")
        ax.set_ylabel("Peaks Reproduced (%)")
        ax.set(ylim=(0, 100))
        ax.xaxis.set_tick_params(labelrotation=45)
        fig.suptitle("Peak Reproducibility")

        return fig, df_data

    # 7d - Fragments within peaks
    def frags_in_peaks(self):
        fig, ax = plt.subplots()

        # Subset data
        df_data = self.metadata_noctrl_table.loc[:, ('id', 'group', 'frip')]

        ax = sns.boxplot(data=df_data, x='group', y='frip', palette = "magma")
        ax.set_ylabel("Fragments within Peaks (%)")
        ax.set(ylim=(0, 100))
        ax.xaxis.set_tick_params(labelrotation=45)
        fig.suptitle("Aligned Fragments within Peaks")

        return fig, df_data
