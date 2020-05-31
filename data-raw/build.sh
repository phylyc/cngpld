#!/bin/sh

make genome=hg19
make genome=hg38
Rscript build.R
