#!/bin/bash

wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/meta_table.csv
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/h3k27me3_R1_raw.csv
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/h3k27me3_R1.bed.stringent.bed
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/h3k27me3_R1.frag_len.txt
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/h3k27me3_R1.frags.bin500.awk.bed
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/h3k27me3_R1.sorted.bam
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/h3k27me3_R2_raw.csv

wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/igg_R1_raw.csv
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/igg_R1.frag_len.txt
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/igg_R1.frags.bin500.awk.bed
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/igg_R1.sorted.bam
wget -nc -P dev/docker/static_reports/test_data/ https://raw.githubusercontent.com/luslab/test-datasets/cutandrun/testdata/dev_reporting/small/igg_R2_raw.csv

docker run --rm -v "$PWD":/home/repo -it luslab/cutandrun-dev-reporting:latest /home/repo/bin/reporting.py gen_reports \
--meta /home/repo/dev/docker/static_reports/test_data/meta_table.csv \
--raw_frag /home/repo/dev/docker/static_reports/test_data/*frag_len.txt \
--bin_frag /home/repo/dev/docker/static_reports/test_data/*bin500*.bed \
--seacr_bed /home/repo/dev/docker/static_reports/test_data/*bed.*.bed \
--bams /home/repo/dev/docker/static_reports/test_data/*sorted.bam \
--output /home/repo/dev/docker/static_reports/test_output \
--log /home/repo/dev/docker/static_reports/test_output/log.txt