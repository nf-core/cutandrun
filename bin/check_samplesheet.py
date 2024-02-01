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
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/cutandrun samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT> <USE_CONTROL>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument(
        "USE_CONTROL",
        help="Boolean for whether or not the user has specified the pipeline must normalise against a control",
    )
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out, use_control):
    """
    This function checks that the samplesheet follows the following structure:

    group,replicate,fastq_1,fastq_2,control
    WT,1,WT_LIB1_REP1_1.fastq.gz,WT_LIB1_REP1_2.fastq.gz,CONTROL_GROUP
    WT,1,WT_LIB2_REP1_1.fastq.gz,WT_LIB2_REP1_2.fastq.gz,CONTROL_GROUP
    WT,2,WT_LIB1_REP2_1.fastq.gz,WT_LIB1_REP2_2.fastq.gz,CONTROL_GROUP
    KO,1,KO_LIB1_REP1_1.fastq.gz,KO_LIB1_REP1_2.fastq.gz,CONTROL_GROUP
    CONTROL_GROUP,1,KO_LIB1_REP1_1.fastq.gz,IGG_LIB1_REP1_2.fastq.gz,
    CONTROL_GROUP,2,KO_LIB1_REP1_1.fastq.gz,IGG_LIB1_REP1_2.fastq.gz,
    """

    # Init
    control_present = False
    num_fastq_list = []
    sample_names_list = []
    control_names_list = []
    sample_run_dict = {}

    with open(file_in, "r") as fin:
        ## Check header
        MIN_COLS = 3
        LEGACY_HEADER = ["group", "replicate", "control_group", "fastq_1", "fastq_2"]
        HEADER = ["group", "replicate", "fastq_1", "fastq_2", "control"]
        HEADER_LEN = len(HEADER)
        header = [x.strip('"') for x in fin.readline().strip().split(",")]

        if len(header) >= len(LEGACY_HEADER) and header[: len(LEGACY_HEADER)] == LEGACY_HEADER:
            print(
                "ERROR: It looks like you are using a legacy header format with a newer version of the pipeline -> {} != {}".format(
                    ",".join(header), ",".join(HEADER)
                )
            )
            sys.exit(1)

        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        line_no = 1
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            ## Check if its just a blank line so we dont error
            if line.strip() == "":
                continue

            ## Check valid number of columns per row
            if len(lspl) != HEADER_LEN:
                print_error(
                    "Invalid number of columns (found {} should be {})! - line no. {}".format(
                        len(lspl), len(HEADER), line_no
                    ),
                    "Line",
                    line,
                )

            ## Set control_present to true if the control column is not empty
            if lspl[4] != "":
                control_present = True

            ## Check valid number of populated columns per row
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, replicate, fastq_1, fastq_2, control = lspl[: len(HEADER)]
            if sample:
                if sample.find(" ") != -1:
                    print_error("Group entry contains spaces!", "Line", line)
            else:
                print_error("Group entry has not been specified!", "Line", line)

            if control:
                if control.find(" ") != -1:
                    print_error("Control entry contains spaces!", "Line", line)

            ## Check for single-end
            if fastq_2 == "":
                print_error("Single-end detected. This pipeline does not support single-end reads!", "Line", line)

            ## Check control sample name is not equal to sample name entry
            if sample == control:
                print_error("Control entry and sample entry must be different!", "Line", line)

            ## Check replicate entry is integer
            if not replicate.isdigit():
                print_error("Replicate id not an integer", "Line", line)
            replicate = int(replicate)
            if replicate <= 0:
                print_error("Replicate must be > 0", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )
            num_fastq = len([fastq for fastq in [fastq_1, fastq_2] if fastq])
            num_fastq_list.append(num_fastq)

            ## Auto-detect paired-end/single-end
            sample_info = []
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = [sample, str(replicate), control, "0", fastq_1, fastq_2]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = [sample, str(replicate), control, "1", fastq_1, fastq_2]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = {sample: {replicate : [ single_end, fastq_1, fastq_2 ]}}
            if sample not in sample_run_dict:
                sample_run_dict[sample] = {}
            if replicate not in sample_run_dict[sample]:
                sample_run_dict[sample][replicate] = [sample_info]
            else:
                if sample_info in sample_run_dict[sample][replicate]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_run_dict[sample][replicate].append(sample_info)

            ## Store unique sample names
            if sample not in sample_names_list:
                sample_names_list.append(sample)

            ## Store unique control names
            if control not in control_names_list:
                control_names_list.append(control)

            line_no = line_no + 1

    ## Check data is either paired-end/single-end and not both
    if min(num_fastq_list) != max(num_fastq_list):
        print_error("Mixture of paired-end and single-end reads!")

    ## Check control group exists
    for ctrl in control_names_list:
        if ctrl != "" and ctrl not in sample_names_list:
            print_error(
                "Each control entry must match at least one group entry! Unmatched control entry: {}.".format(ctrl)
            )

    ## Create control identity variable
    for sample in sorted(sample_run_dict.keys()):
        for replicate in sorted(sample_run_dict[sample].keys()):
            for idx, sample_info in enumerate(sample_run_dict[sample][replicate]):
                if control_present:
                    if sample_info[0] in control_names_list:
                        sample_info.append("1")
                        if sample_info[2] != "":
                            print_error("Control cannot have a control: {}.".format(sample_info[0]))
                    else:
                        sample_info.append("0")
                else:
                    sample_info.append("0")

    ## Check use_control parameter is consistent with input groups
    if use_control == "true" and not control_present:
        print_error(
            "ERROR: No 'control' group was found in "
            + str(file_in)
            + " If you are not supplying a control, please specify --use_control 'false' on command line."
        )

    if use_control == "false" and control_present:
        print(
            "WARNING: Parameter --use_control was set to false, but an control group was found in " + str(file_in) + "."
        )

    # Calculate the exact control/replicate id combo
    if use_control == "true":
        for group, reps in sorted(sample_run_dict.items()):
            # Calculate the ctrl group
            ctrl_group = None
            is_ctrl = False
            for replicate, info in sorted(reps.items()):
                ctrl_group = info[0][2]
                if ctrl_group == "":
                    is_ctrl = True
                break

            # Continue if ctrl
            if is_ctrl:
                continue

            # Get num reps
            num_reps = len(reps)
            num_ctrl_reps = len(sample_run_dict[ctrl_group])

            # Assign actual ctrl rep id
            for rep, info in sorted(reps.items()):
                if num_reps == num_ctrl_reps:
                    ctrl_group_new = ctrl_group + "_" + str(rep)
                else:
                    ctrl_group_new = ctrl_group + "_1"
                info[0][2] = ctrl_group_new

    ## Write validated samplesheet with appropriate columns
    if len(sample_run_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(["id", "group", "replicate", "control", "single_end", "fastq_1", "fastq_2", "is_control"])
                + "\n"
            )
            for sample in sorted(sample_run_dict.keys()):
                ## Check that replicate ids are in format 1..<NUM_REPS>
                uniq_rep_ids = set(sample_run_dict[sample].keys())
                if len(uniq_rep_ids) != max(uniq_rep_ids):
                    print_error(
                        "Replicate ids must start with 1!",
                        "Group",
                        sample,
                    )
                for replicate in sorted(sample_run_dict[sample].keys()):
                    ## Check tech reps have same control group id
                    check_group = sample_run_dict[sample][replicate][0][2]
                    for tech_rep in sample_run_dict[sample][replicate]:
                        if tech_rep[2] != check_group:
                            tech_rep[2] = check_group
                            # print_error("Control group must match within technical replicates", tech_rep[2])

                    ## Write to file
                    for idx, sample_info in enumerate(sample_run_dict[sample][replicate]):
                        sample_id = "{}_R{}_T{}".format(sample, replicate, idx + 1)
                        fout.write(",".join([sample_id] + sample_info) + "\n")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT, args.USE_CONTROL)


if __name__ == "__main__":
    sys.exit(main())
