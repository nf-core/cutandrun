#!/usr/bin/env python

import os
import glob
import argparse

import deeptools.countReadsPerBin as crpb
import pysam

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################
Description = 'Calclate FRIP scores (FRagment proportion in Peaks regions) using deeptools for each sample'

parser = argparse.ArgumentParser(description=Description)

## REQUIRED PARAMETERS
parser.add_argument('--bams', help="Bam file.")
parser.add_argument('--peaks', help="Peaks interval file.")
parser.add_argument('--threads', help="the number of threads for the task.")
parser.add_argument('--outpath', help="Full path to output directory.")
args = parser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

# https://deeptools.readthedocs.io/en/develop/content/example_api_tutorial.html

# Create file lists
bam_file_list = glob.glob(args.bams)
peak_file_list = glob.glob(args.peaks)

frips = []
for idx, bam_file in enumerate(bam_file_list):
    # Init
    frip = 0

    # Read first line
    first_line = None
    with open(peak_file_list[idx], "r") as file:
        for line in file:
            first_line = line
            break

    if first_line is not None:
        print("Calculating " + bam_file + " using " + peak_file_list[idx])
        cr = crpb.CountReadsPerBin([bam_file], bedFile=[peak_file_list[idx]], numberOfProcessors=int(1))

        # Calc the total number of reads in peaks per bam file
        reads_at_peaks = cr.run()
        total = reads_at_peaks.sum(axis=0)

        # Load up bam file and get the total number of mapped reads
        bam = pysam.AlignmentFile(bam_file)

        # Calc frip
        frip = float(total[0]) / bam.mapped

    frips.append(str(frip))

    # Log
    print("Frip = " + str(frip))

# Create string and write to file
frip_string = ",".join(frips)
writer = open(os.path.join(args.outpath, "frips.csv"), "w")
writer.write("frip\n")
writer.write(frip_string)
writer.close()
