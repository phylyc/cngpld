#!/bin/bash

Rscript ../scripts/run.cngpld.R \
  --dir="." \
  --study="luad_vs_lusc"

Rscript ../scripts/plot.cngpld.R \
  --dir="." \
  --study="luad_vs_lusc" \
  --case_label="Lung Adenocarcinoma" \
  --control_label="Lung Squamous Cell Carcinoma" \
  --drivers_dir="."
