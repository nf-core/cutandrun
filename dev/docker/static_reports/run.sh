#!/bin/bash

docker run --rm -v "$PWD":/home/repo -it luslab/cutandrun-dev-reporting:latest /home/repo/bin/python/reporting/main.py genimg \
--input /home/repo/dev/docker/static_reports/test_data/meta_table.csv \
--output /home/repo/dev/docker/static_reports/test_output \
--log /home/repo/dev/docker/static_reports/test_output/log.txt