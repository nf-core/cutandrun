#!/bin/bash

docker run -it -e PASSWORD=password -p 8787:8787 -v ~/dev/repos/cutandrun:/home/rstudio luslab/cutrun-ds2-dev:latest
