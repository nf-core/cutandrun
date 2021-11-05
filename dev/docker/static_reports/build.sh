#!/bin/bash

docker build -f dev/docker/static_reports/Dockerfile -t luslab/cutandrun-dev-reporting:latest dev/docker/static_reports
