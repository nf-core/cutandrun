#!/bin/bash

docker run --rm -v "$PWD":/home/repo -it luslab/cutandrun-dev-reporting:latest /home/repo/bin/python/reporting/main.py genimg \
--meta /home/repo/dev/docker/static_reports/test_data/meta_table.csv \
--raw_frag /home/repo/dev/docker/static_reports/test_data/*_raw.csv \
--output /home/repo/dev/docker/static_reports/test_output \
--log /home/repo/dev/docker/static_reports/test_output/log.txt