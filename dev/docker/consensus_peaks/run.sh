#!/bin/bash

docker run --rm -v "$PWD":/home/repo -it luslab/cutandrun-dev-plot-consensus-peaks:latest /home/repo/bin/consensus_peaks.py \
--peaks /home/repo/dev/docker/consensus_peaks/test_data/*peaks.bed \
--outpath /home/repo/dev/docker/consensus_peaks/test_output
