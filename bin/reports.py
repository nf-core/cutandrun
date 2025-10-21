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

import argparse
import logging
import os
import glob
import pandas as pd

# *
# ========================================================================================
# UTIL
# ========================================================================================
# */


def init_logger(app_name, log_file=None):
    logger = logging.getLogger(app_name)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s:%(name)s:%(levelname)s - %(message)s")

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger


# *
# ========================================================================================
# MAIN
# ========================================================================================
# */


# Input a collection of file tables from different samples and append the sample ids, then output one table
def merge_samples(args):
    # Init
    metadata_path = os.path.abspath(args.metadata)

    # Log
    logger.info("merge_samples " + metadata_path)

    # Get file list
    file_list = glob.glob(metadata_path)
    file_list.sort()

    df_metadata = None
    for idx, file in enumerate(file_list):
        # Strip sample id and group name
        sample_id = os.path.basename(file).replace(args.id_parse_string, "")
        group_name = "_".join(sample_id.split("_")[:-1])

        # Load table
        df_newdata = pd.read_csv(file, sep=",")
        df_newdata["id"] = sample_id
        df_newdata["group"] = group_name
        df_newdata = df_newdata.set_index("id")

        if idx == 0:
            df_metadata = df_newdata
        else:
            df_metadata = df_metadata.append(df_newdata)

    df_metadata.to_csv(args.output, index=True, sep=",")


if __name__ == "__main__":
    # Create command args
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="sub-command help")
    subparsers.required = True

    # Generate reporting function
    newparser = subparsers.add_parser("merge_samples")
    newparser.set_defaults(func=merge_samples)
    newparser.add_argument("--log", required=False)
    newparser.add_argument("--metadata", required=True)
    newparser.add_argument("--id_parse_string", required=True)
    newparser.add_argument("--output", required=True)

    # Parse
    parsed_args = parser.parse_args()

    # Init logging
    logger = init_logger("reporting", parsed_args.log)
    logger.info("CUT&RUN Python Reporting")

    # Call functions
    parsed_args.func(parsed_args)
