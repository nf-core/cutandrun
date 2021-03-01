#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sb
import plotly.express as px

class Figures:
    data_table = None

    def __init__(self, logger, data):
        self.logger = logger
        self.data_path = data

    def load_data(self):
        self.data_table = pd.read_csv(self.data_path, sep=',')


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
        sb.set_theme()
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

    # ---------- Plot 1 - Alignment Summary --------- #
    def alignment_summary(self):
        # Subset data 
        df_data = self.data_table.loc[:, ('id', 'group', 'bt2_total_reads_target', 'bt2_total_aligned_target', 'target_alignment_rate', 'spikein_alignment_rate')]

        ## Construct quad plot
        fig, seq_summary = plt.subplots(2,2)
        fig.suptitle("Sequencing and Alignment Summary")

        # Seq depth
        sb.boxplot(data=df_data, x='group', y='bt2_total_reads_target', order=['h3k27me3', 'h3k4me3', 'igg'], ax=seq_summary[0,0])
        seq_summary[0,0].set_title("Sequencing Depth")
        seq_summary[0,0].set_ylabel("Total Reads")

        # Alignable fragments
        sb.boxplot(data=df_data, x='group', y='bt2_total_aligned_target', order=['h3k27me3', 'h3k4me3', 'igg'], ax=seq_summary[0,1])
        seq_summary[0,1].set_title("Alignable Fragments")
        seq_summary[0,1].set_ylabel("Total Aligned Reads")

        # Alignment rate hg38
        sb.boxplot(data=df_data, x='group', y='target_alignment_rate', order=['h3k27me3', 'h3k4me3', 'igg'], ax=seq_summary[1,0])
        seq_summary[1,0].set_title("Alignment Rate (hg38)")
        seq_summary[1,0].set_ylabel("Percent of Fragments Aligned")

        # Alignment rate e.coli
        sb.boxplot(data=df_data, x='group', y='spikein_alignment_rate', order=['h3k27me3', 'h3k4me3', 'igg'], ax=seq_summary[1,1])
        seq_summary[1,1].set_title("Alignment Rate (e.coli)")
        seq_summary[1,1].set_ylabel("Percent of Fragments Aligned")

        plt.subplots_adjust(wspace=0.5, hspace=0.5)
        
        return fig, df_data

    def alignment_summary_ex(self):
        # Subset data 
        df_data = self.data_table.loc[:, ('id', 'group', 'bt2_total_reads_target', 'bt2_total_aligned_target', 'target_alignment_rate', 'spikein_alignment_rate')]

        fig = px.box(df_data, x="group", y="bt2_total_reads_target")

        return fig, df_data


        
# import glob
# 
# import re


# # ---------- Plots 2&3 - Mapped Fragment Distribution --------- #

# # Parse './' to fragments argument as we expect all data to be in the same directory. 
# fragments = args.exp_fragments

# # Create list of deeptools raw fragment files
# dt_frag_list = glob.glob(fragments + '*raw.csv')

# df_list = list()
# for i in list(range(0,len(dt_frag_list))):
#     # create dataframe from csv file for each file and save to a list
#     dt_frag_i = pd.read_csv(dt_frag_list[i], sep='\t', skiprows=[0], header=0)
#     df_list.append( dt_frag_i )

#     # create long forms of fragment histograms
#     dt_frag_i_long = np.repeat(dt_frag_i['Size'].values, dt_frag_i['Occurrences'].values)
#     dt_sample_i_long = np.repeat(dt_frag_i['Sample'][0], len(dt_frag_i_long))

#     if i==0:
#         frags_arr = dt_frag_i_long
#         sample_arr = dt_sample_i_long
#         og_frag_df = dt_frag_i
#     else:
#         frags_arr = np.append(frags_arr, dt_frag_i_long)
#         sample_arr = np.append(sample_arr, dt_sample_i_long)
#         og_frag_df = og_frag_df.append(dt_frag_i)

# # create hue array using regex pattern matching
# for i in list(range(0,len(sample_arr))):
#     sample_i = sample_arr[i]
#     sample_exp = re.findall("^[^_]*", sample_i)

#     if i==0:
#         sample_exp_arr = np.array(sample_exp[0])
#     else:
#         sample_exp_arr = np.append(sample_exp_arr, sample_exp[0])

# df_long = pd.DataFrame( { "fragment_size" : frags_arr, "sample" : sample_arr , "histone_mark": sample_exp_arr}, index = np.arange(len(frags_arr)))
# df_long.to_csv('./fragmanet_distribution_violin.csv', index=False)
# og_frag_df.to_csv('./fragmanet_distribution_line.csv', index=False)
# # print(df_long)

# plt.clf()
# ax1 = sns.violinplot(data=df_long, x="sample", y="fragment_size", hue="histone_mark")
# # plt.show()
# fig1 = ax1.get_figure()
# fig1.savefig('fragmanet_distribution_violin.png')

# plt.clf()
# ax2 = sns.lineplot(data=og_frag_df, x="Size", y="Occurrences", hue="Sample")
# # plt.show()
# fig2 = ax2.get_figure()
# fig2.savefig('fragmanet_distribution_line.png')

# import datetime
# import numpy as np
# from matplotlib.backends.backend_pdf import PdfPages
# import matplotlib.pyplot as plt

# # Create the PdfPages object to which we will save the pages:
# # The with statement makes sure that the PdfPages object is closed properly at
# # the end of the block, even if an Exception occurs.
# with PdfPages('multipage_pdf.pdf') as pdf:
#     plt.figure(figsize=(3, 3))
#     plt.plot(range(7), [3, 1, 4, 1, 5, 9, 2], 'r-o')
#     plt.title('Page One')
#     pdf.savefig()  # saves the current figure into a pdf page
#     plt.close()

#     # if LaTeX is not installed or error caught, change to `usetex=False`
#     plt.rc('text', usetex=True)
#     plt.figure(figsize=(8, 6))
#     x = np.arange(0, 5, 0.1)
#     plt.plot(x, np.sin(x), 'b-')
#     plt.title('Page Two')
#     pdf.attach_note("plot of sin(x)")  # you can add a pdf note to
#                                        # attach metadata to a page
#     pdf.savefig()
#     plt.close()

#     plt.rc('text', usetex=False)
#     fig = plt.figure(figsize=(4, 5))
#     plt.plot(x, x ** 2, 'ko')
#     plt.title('Page Three')
#     pdf.savefig(fig)  # or you can pass a Figure object to pdf.savefig
#     plt.close()

#     # We can also set the file's metadata via the PdfPages object:
#     d = pdf.infodict()
#     d['Title'] = 'Multipage PDF Example'
#     d['Author'] = 'Jouni K. Sepp\xe4nen'
#     d['Subject'] = 'How to create a multipage pdf file and set its metadata'
#     d['Keywords'] = 'PdfPages multipage keywords author title subject'
#     d['CreationDate'] = datetime.datetime(2009, 11, 13)
#     d['ModDate'] = datetime.datetime.today()