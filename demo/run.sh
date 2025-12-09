#!/bin/bash

Rscript run.cngpld.R \
  --dir="." \
  --study="luad_vs_lusc"

Rscript plot.cngpld.R \
  --dir="." \
  --study="luad_vs_lusc" \
  --case_label="Lung Adenocarcinoma" \
  --control_label="Lung Squamous Cell Carcinoma" \
  --drivers_dir="."
