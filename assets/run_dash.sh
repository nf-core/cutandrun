#!/bin/bash

docker run -it -p 8050:8050 -v "$PWD":/app/data --rm luslab/dash:latest