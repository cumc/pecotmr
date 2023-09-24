# fasterscale

R package for efficiently centering and scaling a matrix.

## Quick Start

To install the latest version of the fasterscale package
from GitHub, use [devtools][devtools]:

```R
install.packages("devtools")
devtools::install_github("pcarbo/fasterscale",build_vignettes = TRUE)
```

This command should automatically install all required packages if
they are not installed already.

If you have cloned the repository locally, you can install the package
with the `install_local` function from devtools. Assuming your working
directory contains the fasterscale repository, run this code to
install the package:

```R
devtools::install_local("fasterscale",build_vignettes = TRUE)
```

[devtools]: https://github.com/r-lib/devtools
