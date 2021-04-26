#!/usr/bin/env python
# coding: utf-8

import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
import seaborn as sns
import plotly.express as px
import pyranges as pr

class Figures:
    data_table = None
    frag_hist = None
    frag_violin = None
    frag_bin500 = None
    seacr_beds = None

    def __init__(self, logger, meta, frags, bin_frag, seacr_bed):
        self.logger = logger
        self.meta_path = meta
        self.frag_path = frags
        self.bin_frag_path = bin_frag
        self.seacr_bed_path = seacr_bed

        sns.set()
        sns.set_theme()
        sns.set_context("paper")

    def format_millions(self, x, pos):
        #the two args are the value and tick position
        return '%1.1fM' % (x * 1e-6)

    def format_thousands(self, x, pos):
        #the two args are the value and tick position
        return '%1.1fK' % (x * 1e-3)

    def load_data(self):
        # ---------- Data - data_table --------- #
        self.data_table = pd.read_csv(self.meta_path, sep=',')

        # ---------- Data - frag_hist --------- #
        # Create list of deeptools raw fragment files
        dt_frag_list = glob.glob(self.frag_path)

        for i in list(range(len(dt_frag_list))):
            # create dataframe from csv file for each file and save to a list
            dt_frag_i = pd.read_csv(dt_frag_list[i], sep='\t', skiprows=[0], header=0)

            # create long forms of fragment histograms
            dt_frag_i_long = np.repeat(dt_frag_i['Size'].values, dt_frag_i['Occurrences'].values)
            dt_sample_i_long = np.repeat(dt_frag_i['Sample'][0], len(dt_frag_i_long))

            if i==0:
                frags_arr = dt_frag_i_long
                sample_arr = dt_sample_i_long
                self.frag_hist = dt_frag_i
            else:
                frags_arr = np.append(frags_arr, dt_frag_i_long)
                sample_arr = np.append(sample_arr, dt_sample_i_long)
                self.frag_hist = self.frag_hist.append(dt_frag_i)

        # ---------- Data - frag_violin --------- #
        # create hue array using regex pattern matching
        for i in list(range(0,len(sample_arr))):
            sample_i = sample_arr[i]
            sample_exp = re.findall("^[^_]*", sample_i)

            if i==0:
                sample_exp_arr = np.array(sample_exp[0])
            else:
                sample_exp_arr = np.append(sample_exp_arr, sample_exp[0])

        self.frag_violin = pd.DataFrame( { "fragment_size" : frags_arr, "sample" : sample_arr , "histone_mark": sample_exp_arr}, index = np.arange(len(frags_arr)))

        # ---------- Data - frag_bin500 --------- #
        # create full join data frame for count data
        # start by creating list of bin500 files
        dt_bin_frag_list = glob.glob(self.bin_frag_path)
        for i in list(range(len(dt_bin_frag_list))):
            dt_bin_frag_i_read = pd.read_csv(dt_bin_frag_list[i], sep='\t', header=None, names=['chrom','bin','count','sample'])
            sample_name = dt_bin_frag_i_read['sample'].iloc[0].split(".")[0]
            dt_bin_frag_i = dt_bin_frag_i_read[['chrom','bin','count']]
            dt_bin_frag_i.columns = ['chrom','bin',sample_name]
            # print(dt_bin_frag_i.head(5))

            if i==0:
                self.frag_bin500 = dt_bin_frag_i
            else:
                self.frag_bin500 = pd.merge(self.frag_bin500, dt_bin_frag_i, on=['chrom','bin'], how='outer')
                # print(self.frag_bin500.head(5))

        # add log2 transformed count data column
        log2_counts = self.frag_bin500[self.frag_bin500.columns[-(len(dt_bin_frag_list)):]].transform(lambda x: np.log2(x))
        chrom_bin_cols = self.frag_bin500[['chrom','bin']]
        self.frag_bin500 = pd.concat([chrom_bin_cols,log2_counts], axis=1)

        # ---------- Data - seacr_bed --------- #
        # create dataframe for seacr peaks
        seacr_bed_list = glob.glob(self.seacr_bed_path)

        # combine all seacr bed files into one df including group and replicate info
        for i in list(range(len(seacr_bed_list))):
            seacr_bed_i = pd.read_csv(seacr_bed_list[i], sep='\t', header=None, usecols=[0,1,2,3,4], names=['chrom','start','end','total_signal','max_signal'])
            bed_base_i = os.path.basename(seacr_bed_list[i])
            sample_id = bed_base_i.split(".")[0]
            [group_i,rep_i] = sample_id.split("_")
            seacr_bed_i['group'] = np.repeat(group_i, seacr_bed_i.shape[0])
            seacr_bed_i['replicate'] = np.repeat(rep_i, seacr_bed_i.shape[0])

            if i==0:
                self.seacr_beds = seacr_bed_i

            else:
                self.seacr_beds = self.seacr_beds.append(seacr_bed_i)

    def annotate_data(self):
        # Make new perctenage alignment columns
        self.data_table['target_alignment_rate'] = self.data_table.loc[:, ('bt2_total_aligned_target')] / self.data_table.loc[:, ('bt2_total_reads_target')] * 100
        self.data_table['spikein_alignment_rate'] = self.data_table.loc[:, ('bt2_total_aligned_spikein')] / self.data_table.loc[:, ('bt2_total_reads_spikein')] * 100
        # self.data_table.describe()
        # print(self.data_table)
        # self.data_table.info()

    def generate_plots(self):
        # Init
        plots = dict()
        data = dict()

        # Get Data
        self.load_data()
        self.annotate_data()

        # Plot 1
        plot1, data1 = self.alignment_summary()
        plots["alignment_summary"] = plot1
        data["alignment_summary"] = data1

        # Plot 2
        plot2, data2 = self.duplication_summary()
        plots["duplication_summary"] = plot2
        data["duplication_summary"] = data2

        # Plot 3
        plot3, data3 = self.fraglen_summary_violin()
        plots["frag_violin"] = plot3
        data["frag_violin"] = data3

        # Plot 4
        plot4, data4 = self.fraglen_summary_histogram()
        plots["frag_hist"] = plot4
        data["frag_hist"] = data4

        # Plot 5
        plot5, data5 = self.replicate_heatmap()
        plots["replicate_heatmap"] = plot5
        data["replicate_heatmap"] = data5

        # Plot 6
        plot6, data6 = self.scale_factor_summary()
        plots["scale_factor_summary"] = plot6
        data["scale_factor_summary"] = data6
        
        # Plot 7a 
        plot7a, data7a = self.no_of_peaks()
        plots["no_of_peaks"] = plot7a
        data["no_of_peaks"] = data7a

        # Plot 7b
        plot7b, data7b = self.peak_widths()
        plots["peak_widths"] = plot7b
        data["peak_widths"] = data7b

        # Plot 7c 
        plot7c, data7c = self.reproduced_peaks()
        plots["reproduced_peaks"] = plot7c
        data["reproduced_peaks"] = data7c

        return (plots, data)

    def generate_dash_plots(self):
        # Init
        plots = dict()
        data = dict()

        # Get Data
        self.load_data()
        self.annotate_data()

        # Plot 1
        plot1, data1 = self.alignment_summary_ex()
        plots["alignment_summary"] = plot1
        data["alignment_summary"] = data1

        return (plots, data)

    def gen_plots_to_folder(self, output_path):
        # Init
        abs_path = os.path.abspath(output_path)

        # Get plots and supporting data tables
        plots, data = self.generate_plots()

        # Save data to output folder
        for key in data:
            data[key].to_csv(os.path.join(abs_path, key + '.csv'), index=False)
            plots[key].savefig(os.path.join(abs_path, key + '.png'))

        # Save pdf of the plots
        self.gen_pdf(abs_path, plots)

    def gen_pdf(self, output_path, plots):
        with PdfPages(os.path.join(output_path, 'report.pdf')) as pdf:
            for key in plots:
                pdf.savefig(plots[key])

            # # We can also set the file's metadata via the PdfPages object:
            # d = pdf.infodict()
            # d['Title'] = 'Multipage PDF Example'
            # d['Author'] = 'Jouni K. Sepp\xe4nen'
            #        d['Subject'] = 'How to create a multipage pdf file and set its metadata'
            # d['Keywords'] = 'PdfPages multipage keywords author title subject'
            # d['CreationDate'] = datetime.datetime(2009, 11, 13)
            # d['ModDate'] = datetime.datetime.today()

    ##### PLOTS #####

    # TODO: Add group ordering -  order=['h3k27me3', 'h3k4me3', 'igg']

    # ---------- Plot 1 - Alignment Summary --------- #
    def alignment_summary(self):
        # Subset data 
        df_data = self.data_table.loc[:, ('id', 'group', 'bt2_total_reads_target', 'bt2_total_aligned_target', 'target_alignment_rate', 'spikein_alignment_rate')]

        ## Construct quad plot
        fig, seq_summary = plt.subplots(2,2)
        fig.suptitle("Sequencing and Alignment Summary")

        # Seq depth
        sns.boxplot(data=df_data, x='group', y='bt2_total_reads_target', ax=seq_summary[0,0])
        seq_summary[0,0].set_title("Sequencing Depth")
        seq_summary[0,0].set_ylabel("Total Reads")

        # Alignable fragments
        sns.boxplot(data=df_data, x='group', y='bt2_total_aligned_target', ax=seq_summary[0,1])
        seq_summary[0,1].set_title("Alignable Fragments")
        seq_summary[0,1].set_ylabel("Total Aligned Reads")

        # Alignment rate hg38
        sns.boxplot(data=df_data, x='group', y='target_alignment_rate', ax=seq_summary[1,0])
        seq_summary[1,0].set_title("Alignment Rate (hg38)")
        seq_summary[1,0].set_ylabel("Percent of Fragments Aligned")

        # Alignment rate e.coli
        sns.boxplot(data=df_data, x='group', y='spikein_alignment_rate', ax=seq_summary[1,1])
        seq_summary[1,1].set_title("Alignment Rate (e.coli)")
        seq_summary[1,1].set_ylabel("Percent of Fragments Aligned")

        plt.subplots_adjust(wspace=0.5, hspace=0.5)
        
        return fig, df_data

    # ---------- Plot 2 - Duplication Summary --------- #
    def duplication_summary(self):
        # Init
        k_formatter = FuncFormatter(self.format_thousands)
        m_formatter = FuncFormatter(self.format_millions)

        # Subset data 
        df_data = self.data_table.loc[:, ('id', 'group', 'dedup_percent_duplication', 'dedup_estimated_library_size', 'dedup_read_pairs_examined')]
        df_data['dedup_percent_duplication'] *= 100
        df_data['unique_frag_num'] = df_data['dedup_read_pairs_examined'] * (1-df_data['dedup_percent_duplication']/100)

        ## Construct quad plot
        fig, seq_summary = plt.subplots(1,3)
        fig.suptitle("Duplication Summary")

        # Duplication rate
        sns.boxplot(data=df_data, x='group', y='dedup_percent_duplication', ax=seq_summary[0])
        seq_summary[0].set_ylabel("Duplication Rate (%)")
        seq_summary[0].set(ylim=(0, 100))
        seq_summary[0].xaxis.set_tick_params(labelrotation=45)

        # Estimated library size
        sns.boxplot(data=df_data, x='group', y='dedup_estimated_library_size', ax=seq_summary[1])
        seq_summary[1].set_ylabel("Estimated Library Size")
        seq_summary[1].yaxis.set_major_formatter(m_formatter)
        seq_summary[1].xaxis.set_tick_params(labelrotation=45)

        # No. of unique fragments 
        sns.boxplot(data=df_data, x='group', y='unique_frag_num', ax=seq_summary[2])
        seq_summary[2].set_ylabel("No. of Unique Fragments")
        seq_summary[2].yaxis.set_major_formatter(k_formatter)
        seq_summary[2].xaxis.set_tick_params(labelrotation=45)

        # Set the plot sizing
        plt.subplots_adjust(top = 0.9, bottom = 0.2, right = 0.9, left = 0.1, hspace = 0.7, wspace = 0.7)
        
        return fig, df_data

    
    # ---------- Plot 3 - Fragment Distribution Violin --------- #
    def fraglen_summary_violin(self):
        fig, ax = plt.subplots()
        ax = sns.violinplot(data=self.frag_violin, x="sample", y="fragment_size", hue="histone_mark")

        return fig, self.frag_violin

    # ---------- Plot 4 - Fragment Distribution Histogram --------- #
    def fraglen_summary_histogram(self):
        fig, ax = plt.subplots()
        ax = sns.lineplot(data=self.frag_hist, x="Size", y="Occurrences", hue="Sample")

        return fig, self.frag_hist

    def alignment_summary_ex(self):
        df_data = self.data_table.loc[:, ('id', 'group', 'bt2_total_reads_target', 'bt2_total_aligned_target', 'target_alignment_rate', 'spikein_alignment_rate')]

        ax = px.box(df_data, x="group", y="bt2_total_reads_target")

        return ax, df_data


    # ---------- Plot 5 - Replicate Reproducibility Heatmap --------- #
    def replicate_heatmap(self):
        sns.set(font_scale=0.6)
        fig, ax = plt.subplots()
        plot_data = self.frag_bin500[self.frag_bin500.columns[-(len(self.frag_bin500.columns)-2):]]
        corr_mat = plot_data.corr()
        ax = sns.heatmap(corr_mat, annot=True)

        return fig, self.frag_bin500

    # ---------- Plot 6 - Replicate Reproducibility Heatmap --------- #
    def scale_factor_summary(self):
        fig, scale_summary = plt.subplots(1,2)
        fig.suptitle("Scaling Factor")

        # Get normalised count data
        df_normalised_frags = self.data_table.loc[:, ('id', 'group')]
        df_normalised_frags['normalised_frags'] = self.data_table['bt2_total_reads_target']*self.data_table['scale_factor']

        # Subset meta data
        df_data_scale = self.data_table.loc[:, ('id', 'group','scale_factor')]

        # Scale factor
        sns.boxplot(data=df_data_scale, x='group', y='scale_factor', ax=scale_summary[0])
        scale_summary[0].set_ylabel('Scale Factor')

        # Normalised fragment count
        sns.boxplot(data=df_normalised_frags, x='group', y='normalised_frags', ax=scale_summary[1])
        scale_summary[1].set_ylabel('Normalised Fragment Count')

        return fig, df_data_scale
    
    # ---------- Plot 7 - Peak Analysis --------- #
    def no_of_peaks(self):
        # 7a - Number of peaks
        fig, ax = plt.subplots()
        fig.suptitle("Total Peaks")

        ## create number of peaks df
        unique_groups = self.seacr_beds.group.unique()
        unique_replicates = self.seacr_beds.replicate.unique()
        self.df_no_peaks = pd.DataFrame(index=range(0,(len(unique_groups)*len(unique_replicates))), columns=['group','replicate','all_peaks'])
        k=0 # counter

        for i in list(range(len(unique_groups))):
            for j in list(range(len(unique_replicates))):
                self.df_no_peaks.at[k,'all_peaks'] = self.seacr_beds[(self.seacr_beds['group']==unique_groups[i]) & (self.seacr_beds['replicate']==unique_replicates[j])].shape[0]
                self.df_no_peaks.at[k,'group'] = unique_groups[i]
                self.df_no_peaks.at[k,'replicate'] = unique_replicates[j]
                k=k+1

        ax = sns.boxplot(data=self.df_no_peaks, x='group', y='all_peaks')
        ax.set_ylabel("No. of Peaks")

        return fig, self.df_no_peaks

        # 7b - Width of peaks
    def peak_widths(self):
        fig, ax = plt.subplots()

        ## add peak width column
        self.seacr_beds['peak_width'] = self.seacr_beds['end'] - self.seacr_beds['start']
        self.seacr_beds['peak_width'] = self.seacr_beds['peak_width'].abs()

        ax = sns.violinplot(data=self.seacr_beds, x="group", y="peak_width", hue="replicate")

        return fig, self.seacr_beds


        # 7c - Peaks reproduced
    def reproduced_peaks(self):
        fig, ax = plt.subplots()
        ## empty dataframe to fill in loop
        reprod_peak_stats = self.df_no_peaks
        reprod_peak_stats = reprod_peak_stats.reindex(columns=reprod_peak_stats.columns.tolist() + ['no_peaks_reproduced','peak_reproduced_rate'])

        ## create pyranges objects and fill df
        unique_groups = self.seacr_beds.group.unique()
        unique_replicates = self.seacr_beds.replicate.unique()
        # pyr_overlap = pr.PyRanges(chromosomes="chr1", starts=(1, 5), ends=[3, 149],
        #     strands=("+", "-"), int64=True)
        # print(len(pyr_overlap))

        # peak_ranges = pr.PyRanges(chromosomes=self.seacr_beds['chrom'], starts=self.seacr_beds['start'], ends=self.seacr_beds['end'])
        # print(peak_ranges)
        # fill_no_peaks = np.empty((0,))

        for i in list(range(len(unique_groups))):
            pyr_query = pr.PyRanges()
            group_i = unique_groups[i]
            for j in list(range(len(unique_replicates))):
                rep_i = unique_replicates[j]
                peaks_i = self.seacr_beds[(self.seacr_beds['group']==group_i) & (self.seacr_beds['replicate']==rep_i)]
                pyr_subject = pr.PyRanges(chromosomes=peaks_i['chrom'], starts=peaks_i['start'], ends=peaks_i['end'])
                if(len(pyr_query) > 0):
                    pyr_overlap = pyr_query.join(pyr_subject)
                    pyr_overlap = pyr_overlap.apply(lambda df: df.drop(['Start_b','End_b'], axis=1))
                    pyr_query = pyr_overlap
                else:
                    pyr_query = pyr_subject
            # fill_no_peaks
            
            # idx = ((j)*len(unique_groups)) - (len(unique_groups) - (i+1))
            # print(idx)        
            # reprod_peak_stats.at[idx, 'no_peaks_reproduced'] = len(pyr_overlap)
            # reprod_peak_stats.at[idx+1, 'no_peaks_reproduced'] = len(pyr_overlap)



        return fig, reprod_peak_stats
        # 7d - Fragments within peaks
        
