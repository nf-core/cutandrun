#!/usr/bin/env python
# coding: utf-8

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys
import os
import argparse
import logging

from lib.figures import Figures

def init_logger(app_name, log_file = None):
    logger = logging.getLogger(app_name)
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s:%(name)s:%(levelname)s - %(message)s')

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger

def gen_png(parsed_args):
    meta_path = parsed_args.meta
    frag_path = parsed_args.raw_frag
    output_path = parsed_args.output
    logger = init_logger('gen_img', parsed_args.log)
    bin_frag_path = parsed_args.bed_frag
    
    logger.info('Generating plots to output folder')
    fig = Figures(logger, meta_path, frag_path)
    fig.gen_plots_to_folder(output_path)

    logger.info('Completed')

if __name__ == '__main__': 
    # Create command args
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    # Functions
    parser_genimg = subparsers.add_parser('genimg', help='generate images help')
    parser_genimg.set_defaults(func=gen_png)
    parser_genimg.add_argument('--meta', required=True)
    parser_genimg.add_argument('--raw_frag', required=True)
    parser_genimg.add_argument('--output', required=True)
    parser_genimg.add_argument('--log', required=False)
    parser_genimg.add_argument('--bed_frag', required=True)

    # Parse
    parsed_args = parser.parse_args()

    # Init logging
    logger = init_logger('reporting', parsed_args.log)
    logger.info("Reporting app")

    # Call functions
    parsed_args.func(parsed_args)