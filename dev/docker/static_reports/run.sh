#!/bin/bash

docker run --rm -v "$PWD":/home/repo -it luslab/cutandrun-dev-reporting:latest /home/repo/bin/reporting.py gen_reports \
--meta /home/repo/results/04_reporting/meta_table.csv \
--raw_frag /home/repo/results/03_peak_calling/06_fragments/*frag_len.txt \
--bin_frag /home/repo/results/03_peak_calling/06_fragments/*bin500*.bed \
--seacr_bed /home/repo/results/03_peak_calling/04_called_peaks/*bed.*.bed \
--bams /home/repo/results/02_alignment/bowtie2/target/markdup/*.bam \
--output /home/repo/dev/docker/static_reports/test_output \
--log /home/repo/dev/docker/static_reports/test_output/log.txt
