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

def gen_png(parsed_args):
    input_path = parsed_args.input
    output_path = parsed_args.output
    
    log_path = None
    if parsed_args.log:
        log_path = parsed_args.log
        print(log_path)

    logging.basicConfig(filename=log_path, level=os.environ.get("LOGLEVEL", "INFO"))
    logger = logging.getLogger("gen_png")
    logger.info('Generating plots to output folder')

    fig = Figures(logger, input_path)
    fig.gen_plots_to_folder(output_path)


if __name__ == '__main__':
    logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
    logger = logging.getLogger("app")
    logger.info(sys.executable)
 
    # Create command args
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    # Functions
    parser_genimg = subparsers.add_parser('genimg', help='generate images help')
    parser_genimg.set_defaults(func=gen_png)
    parser_genimg.add_argument('--input', required=True)
    parser_genimg.add_argument('--output', required=True)
    parser_genimg.add_argument('--log', required=False)

    # Parse
    parsed_args = parser.parse_args()
    parsed_args.func(parsed_args)