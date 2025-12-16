#!/bin/bash

Rscript -e 'library(devtools); document(); install();'
Rscript -e 'install.packages("optparse", repos="https://cran.r-project.org")'
