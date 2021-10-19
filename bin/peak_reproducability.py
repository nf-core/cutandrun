#!/usr/bin/env python

import os
import glob
import argparse

import dask.dataframe as dd
import numpy as np
import pandas as pd

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################
Description = 'Calclate peak reproducability percentage for each sample'

parser = argparse.ArgumentParser(description=Description)

## REQUIRED PARAMETERS
parser.add_argument('--intersect', help="Peaks intersect file.")
parser.add_argument('--threads', help="the number of threads for the task.")
parser.add_argument('--outpath', help="Full path to output directory.")
args = parser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

# Init
peak_perc = 0

print('Reading file')

# Read file in using dask
ddf_inter = dd.read_csv(args.intersect, sep='\t', header=None, names=['chrom','start','end','overlap_1','key','a_name','b_name','count'])

# Find number of files
numfiles = ddf_inter['b_name'].max().compute()

# Check for table format
if isinstance(numfiles, str):
    print('Detected single file, reloading table')
    numfiles = 1
    ddf_inter = dd.read_csv(args.intersect, sep='\t', header=None, names=['chrom','start','end','overlap_1','overlap_2','key','name','count'])
    
print('Number of files: ' + str(numfiles))

# Check for empty file
if numfiles != 0:
    # Find total number of peaks
    ddf_inter_grouped = ddf_inter.groupby(by=["key"]).size()
    df_inter_grouped = ddf_inter_grouped.compute()
    total_peaks = len(df_inter_grouped.index)
    print('Total peaks: ' + str(total_peaks))

    if total_peaks > 0:
        # Filter for files which had an overlap and group by peak
        ddf_inter_filt = ddf_inter[ddf_inter["count"] > 0]
        ddf_inter_grouped = ddf_inter_filt.groupby(by=["key"]).size()
        df_inter_grouped = ddf_inter_grouped.compute()
        df_inter_grouped = df_inter_grouped.reset_index()
        df_inter_grouped = df_inter_grouped.rename({0: 'count'}, axis=1)

        # Filter for peaks which have full overlap
        df_inter_grouped_filter = df_inter_grouped[df_inter_grouped["count"] == numfiles]
        overlap_peaks = len(df_inter_grouped_filter.index)
        print('Overlap peaks: ' + str(overlap_peaks))

        # Calc peak percentage
        peak_perc = (overlap_peaks / total_peaks) * 100

# Create string and write to file
output_string = str(peak_perc)
writer = open(os.path.join(args.outpath, "peak_repro.csv"), "w")
writer.write("peak_repro\n")
writer.write(output_string)
writer.close()