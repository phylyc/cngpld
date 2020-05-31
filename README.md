# Case-control copy-number analysis 

Package for performing a case-control copy-number analysis
using Gaussian process latent difference from copy-number seg
files or GISTIC2 outputs. This statistical model identities regions that 
are more frequently disrupted in the case cohort compared to the control.

### Installation ###

Install the [gpldiff](https://bitbucket.org/djhshih/gpldiff) package.

Generated the required data for the package by

```
cd data-raw
./build.sh
```

Then, this package may be installed using `devtools`.
The documentation needs to be generated prior to installation.

```
library(devtools)
document()
install()
```

### Example ###

From within the `demo` directory.

Retrieve and preprocess the input files by
```
make
```

Then, perform the analysis by executing `run-partial.R` or `run.R`.

