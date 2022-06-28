# Case-control copy-number analysis 

Package for performing a case-control copy-number analysis
using Gaussian process latent difference from copy-number seg
files or GISTIC2 outputs. This statistical model identities regions that 
are more frequently disrupted in the case cohort compared to the control.

### Installation ###

1. Install the [gpldiff](https://bitbucket.org/djhshih/gpldiff) package.

2. Install the [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) package from Bioconductor.

3. Clone this repository by `git clone https://bitbucket.org/djhshih/cngpld`

4. Navigate to the directory of this repository in your shell environment. Generate the required data files for the package by

```
cd data-raw
./build.sh
```

You can edit `build.sh` to generate additional coordinates for additional genomes. See `Makefile`.

5. Navigate to the directory of this repository in your R environment. Generate the documentation for required for this package and install the package using `devtools`.

```
library(devtools)
document()
install()
```

Note: Do *not* simply run `devtools::install_bitbucket()` as it will not generated the required data and documentation files.


### Example ###

From within the `demo` directory.

Retrieve and preprocess the input files by
```
make
```

Then, perform the analysis by executing `run-partial.R` or `run.R`.


### Citation ####

If you use this tool in a scientific work , please cite:

CNGPLD: case–control copy-number analysis using Gaussian process latent difference.
David J H Shih, Ruoxing Li, Peter Müller, W Jim Zheng, Kim-Anh Do, Shiaw-Yih Lin, Scott L Carter.
Bioinformatics. 2022 Feb 17;38(8):2096-2101. doi: 10.1093/bioinformatics/btac096.