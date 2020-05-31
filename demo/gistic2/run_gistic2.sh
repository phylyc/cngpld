#!/bin/sh

# Run GISTIC2 on on input seg file

# Prequisite: gistic2 >= 2.0.23 is available on $PATH
# gistic2 may be installed with (patefiant)[https://github.com/djhshih/patefiant]

segfile=$1
outdir=$2

gistic_path=$(which gistic2)
gistic_root=${gistic_path%/bin/gistic2}

mkdir -p $outdir

# input files
refgenefile=$gistic_root/refgenefiles/hg19.mat
# markers file is not required in v2.0.23
# cnv file will be ignored without the markers file

gistic2 -b $outdir -seg $segfile -refgene $refgenefile -cnv $cnvfile -genegistic 1 -smallmem 1 -broad 1 -brlen 0.5 -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme

