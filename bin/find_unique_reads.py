#!/usr/bin/env python

import os
import argparse
import csv

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################
Description = "Find names of the unique alignments"

parser = argparse.ArgumentParser(description=Description)

## REQUIRED PARAMETERS
parser.add_argument("--bed_path")
parser.add_argument("--output_path")
parser.add_argument("--metrics_path")
args = parser.parse_args()

############################################
############################################
## MAIN FUNCTION
############################################
############################################

# Init
output_path = os.path.abspath(args.output_path)
bed_path = os.path.abspath(args.bed_path)
metrics_path = os.path.abspath(args.metrics_path)

# Empty dictionary for the alignments
alignments = dict()

# Counter for alignments before filtering
i = 0

# Open BED file containing new line for each alignment
with open(bed_path) as csvfile:

    # Loop through CSV reader object
    data = csv.reader(csvfile, delimiter="\t")
    for row in data:

        # Skip read 2
        if row[3].endswith("2"):
            continue

        # String for unique reads, starts with chromosome
        unique_id = row[0]

        # If read 1 maps to + strand, use start column as starting point
        if row[5] == "+":
            unique_id += "+" + row[1]
        # If read 2 maps to - strand, use end column as starting point
        else:
            unique_id += "-" + row[2]

        # If the ID is not already in the results, add it
        if unique_id not in alignments:
            alignments[unique_id] = [row[3], row[4]]

        # If ID already exists, check whether mapping quality is better
        elif alignments[unique_id][1] > row[4]:
            alignments[unique_id] = [row[3], row[4]]

        # Increase counter
        i += 1

# Extract alignment names only (QNAME from BAM)
alignments = [i[1][0][:-2] for i in alignments.items()]

# Collect metrics into a string
report = "LINEAR AMPLIFICATION DUPLICATION METRICS"
report += f"\nReads before filtering\t{i}"
report += f"\nLA duplicates removed (n)\t{i-len(alignments)}"
report += f"\nLA duplicates removed (%)\t{round((i-len(alignments))/i*100, 2)}"
report += f"\nUnique reads after LA duplicate removal\t{len(alignments)}"

# Write string to a text file
with open(metrics_path, "w") as f:
    f.write(report)

# Write alignment names to a text file
with open(output_path, "w") as f:
    for alignment in alignments:
        f.write("%s\n" % alignment)
