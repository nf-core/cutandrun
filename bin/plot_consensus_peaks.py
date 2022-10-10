#!/usr/bin/env python

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import upsetplot

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################
Description = "Upset venn diagram of consensus peaks."
Epilog = """Example usage: python consensus_peaks.py <MERGED_INTERVAL_FILE> """

parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
parser.add_argument("--peaks", help="Merged peaks interval file with replicate counts column.")
parser.add_argument("--outpath", help="Full path to output directory.")
args = parser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

# create list of data frames, one for each group consensus peaks file
peak_file_list = glob.glob(args.peaks)

if len(peak_file_list) > 10:
    print("WARN: There are too many files to generate an upset plot, cancelling figure generation")
    exit(0)

peak_df_list = list()
for i in list(range(len(peak_file_list))):
    peaks_i = pd.read_csv(
        peak_file_list[i],
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 8, 9],
        names=["chrom", "start", "end", "sample_reps", "count"],
    )
    peaks_i["sample_reps"] = peaks_i["sample_reps"].replace(".peaks.bed.stringent.bed", "", regex=True)
    peak_df_list.append(peaks_i)
    reps2 = peaks_i[peaks_i["count"] > 1]

# add sorted column to each dataframe, and make new condensed dataframe list
summary_peak_df_list = list()
for i in list(range(len(peak_df_list))):
    peaks_i = peak_df_list[i]
    peaks_i["sorted_samples"] = ""
    rows_now = peaks_i.shape[0]
    for j in list(range(rows_now)):
        sample_list = peaks_i.at[j, "sample_reps"]
        sample_array = np.unique(sample_list.split(","))
        sample_sorted = sorted(sample_array)
        sample_str = ",".join(sample_sorted)
        peaks_i.at[j, "sorted_samples"] = sample_str
    summary_peaks_i = peaks_i[["sorted_samples", "count"]].groupby(["sorted_samples"], as_index=False).sum()
    summary_peak_df_list.append(summary_peaks_i)

# construct data in appropriate format for upsetplot, and plot
for i in list(range(len(summary_peak_df_list))):
    df_i = summary_peak_df_list[i]
    # Get group name
    basename = os.path.basename(peak_file_list[i])
    group_name = basename.rsplit(".", -1)[0]
    file_name = group_name + ".consensus_peaks.pdf"
    categories = df_i.shape[0]
    cat_list = []
    for j in list(range(categories)):
        summary_sample = df_i.at[j, "sorted_samples"].split(",")
        cat_list.append(summary_sample)

    # Plot
    peak_counts = upsetplot.from_memberships(cat_list, data=df_i["count"])
    upsetplot.plot(peak_counts)
    plt.show()
    plt.savefig(os.path.join(args.outpath, file_name))
