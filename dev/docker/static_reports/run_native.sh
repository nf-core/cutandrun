#!/bin/bash

bin/reporting.py gen_reports \
--meta "results/04_reporting/meta_table.csv" \
--raw_frag "results/03_peak_calling/06_fragments/*frag_len.txt" \
--bin_frag "results/03_peak_calling/06_fragments/*bin500*.bed" \
--seacr_bed "results/03_peak_calling/04_called_peaks/*bed.*.bed" \
--bams "results/02_alignment/bowtie2/target/markdup/*.bam" \
--output "dev/docker/static_reports/test_output" \
--log "dev/docker/static_reports/test_output/log.txt"
