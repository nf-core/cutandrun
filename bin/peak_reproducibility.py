#!/usr/bin/env python

# MIT License
#
# Copyright (c) 2023 @chris-cheshire
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Author: @chris-cheshire

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
Description = "Calclate peak reproducibility percentage for each sample"

parser = argparse.ArgumentParser(description=Description)

## REQUIRED PARAMETERS
parser.add_argument("--sample_id", help="Sample id.")
parser.add_argument("--intersect", help="Peaks intersect file.")
parser.add_argument("--threads", help="the number of threads for the task.")
parser.add_argument("--outpath", help="Full path to output directory.")
args = parser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

# Init
peak_perc = 0
numfiles = 0
num_columns = 0

print("Reading file")

# Read first line
first_line = None
with open(args.intersect, "r") as file:
    for line in file:
        first_line = line
        break

if first_line is not None:
    first_line_split = first_line.split("\t")
    num_columns = len(first_line_split)
    numfiles = 1

if numfiles != 0:
    print("Number of columns: " + str(num_columns))

    ddf_inter = None
    if num_columns == 6:
        # Read file in using dask
        ddf_inter = dd.read_csv(
            args.intersect,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "key", "file_num", "count"],
            dtype={
                "chrom": str,
                "start": np.int64,
                "end": np.int64,
                "key": str,
                "file_num": np.int32,
                "count": np.int32,
            },
        )
        numfiles = ddf_inter["file_num"].max().compute()

    elif num_columns == 5:
        # Read file in using dask
        ddf_inter = dd.read_csv(
            args.intersect,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "key", "count"],
            dtype={
                "chrom": str,
                "start": np.int64,
                "end": np.int64,
                "key": str,
                "file_num": np.int32,
                "count": np.int32,
            },
        )
        numfiles = 1
    else:
        print("Invalid file format detected")
        exit(1)

    print("Number of files: " + str(numfiles))

    # Find total number of peaks
    ddf_inter_grouped = ddf_inter.groupby(by=["key"]).size()
    df_inter_grouped = ddf_inter_grouped.compute()
    total_peaks = len(df_inter_grouped.index)
    print("Total peaks: " + str(total_peaks))

    if total_peaks > 0:
        # Filter for files which had an overlap and group by peak
        ddf_inter_filt = ddf_inter[ddf_inter["count"] > 0]
        ddf_inter_grouped = ddf_inter_filt.groupby(by=["key"]).size()
        df_inter_grouped = ddf_inter_grouped.compute()
        df_inter_grouped = df_inter_grouped.reset_index()
        df_inter_grouped = df_inter_grouped.rename({0: "count"}, axis=1)

        # Filter for peaks which have full overlap
        df_inter_grouped_filter = df_inter_grouped[df_inter_grouped["count"] == numfiles]
        overlap_peaks = len(df_inter_grouped_filter.index)
        print("Overlap peaks: " + str(overlap_peaks))

        # Calc peak percentage
        peak_perc = (overlap_peaks / total_peaks) * 100
else:
    print("Empty file detected")

# Create string and write to file
output_string = str(peak_perc)
writer = open(os.path.join(args.outpath, args.sample_id + "_peak_repro.tsv"), "w")
writer.write(args.sample_id + "\t" + output_string + "\n")
writer.close()
