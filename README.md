# pecotmr [![codecov](https://codecov.io/gh/cumc/pecotmr/graph/badge.svg?token=NFH9R1OPG7)](https://codecov.io/gh/cumc/pecotmr)

This R package, `pecotmr` implements **p**air-wise **e**nrichment, **co**localization, **T**WAS and **M**endelian **R**andomization to integrate QTL and GWAS analysis based on fine-mapped single effects model.
`pecotmr` also stands for a **p**robabilistic **eco**system consolidating a range of established and new models into a streamlined, user-friendly interface complete with well-documented examples that will help getting your analysis done by **t**o**m**o**r**row if you start working on it today!

This project was developed by NIH/NIA FunGen-xQTL Project Working Group members with applications to FunGen-xQTL data analysis.

## Quick Start

To install the latest version of the pecotmr package
from GitHub, use [devtools][devtools]:

```R
install.packages("devtools")
devtools::install_github("cumc/pecotmr",build_vignettes = TRUE)
```

This command should automatically install all required packages if
they are not installed already.

If you have cloned the repository locally, you can install the package
with the `install_local` function from devtools. Assuming your working
directory contains the `pecotmr` repository, run this code to
install the package:

```R
devtools::install_local("pecotmr",build_vignettes = TRUE)
```

## Developer's notes

+ When any changes are made to `roxygen2` markup or the C++ code in
the src directory, run `devtools::document()` to update the
[RcppExports.cpp](src/RcppExports.cpp), the package namespaces (see
[NAMESPACE](NAMESPACE)), and the package documentation files (in the
"man" subdirectory),

+ These are the R commands to build the website (make sure you are
connected to Internet while running these commands):

   ```R
   pkgdown::build_site(lazy=TRUE, examples=FALSE)
   ```

+ After editing C++ code in the `src` directory, please use
[uncrustify][uncrustify] (version >=0.74.0, available from conda-forge) 
to format the code using configuration file
`inst/misc/uncrustify_default.cfg`. For example:

   ```bash
   uncrustify -c inst/misc/uncrustify_default.cfg --replace --no-backup -l CPP src/qtl_enrichment.cpp
   uncrustify -c inst/misc/uncrustify_default.cfg --replace --no-backup -l CPP src/qtl_enrichment.hpp
   ```

+ Prior to submitting the package to CRAN, the following modifications
need to be made: (1) remove the `Remotes:` entry in `DESCRIPTION`; and
(2) remove the `fastenloc.Rmd` vignette.

[devtools]: https://github.com/r-lib/devtools
[uncrustify]: https://github.com/uncrustify/uncrustify
