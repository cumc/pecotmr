# pecotmr

[![CI](https://github.com/StatFunGen/pecotmr/actions/workflows/ci.yml/badge.svg?branch=)](https://github.com/StatFunGen/pecotmr/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/StatFunGen/pecotmr/graph/badge.svg?token=NFH9R1OPG7)](https://codecov.io/gh/StatFunGen/pecotmr)
[![Website](https://github.com/StatFunGen/pecotmr/actions/workflows/pkgdown.yml/badge.svg)](https://statfungen.github.io/pecotmr)

This R package, `pecotmr` implements **p**ost-GWAS **e**nrichment, **co**localization, **T**WAS and **M**endelian **R**andomization to integrate QTL and GWAS analysis, mostly based on fine-mapped single effects model although a number of other approaches are included for handy availiability in daily research.
`pecotmr` also stands for a **p**robabilistic **eco**system consolidating a range of established and new models into a streamlined, user-friendly interface complete with well-documented examples that will help getting your analysis done by **t**o**m**o**r**row if you start working on it today!

This project was developed by NIH/NIA FunGen-xQTL Project Working Group members with applications to FunGen-xQTL data analysis.

## Quick Start

To install the latest version of the pecotmr package
from GitHub, use [BiocManager][BiocManager]:

```R
install.packages(c("BiocManager", "remotes"))
BiocManager::install("StatFunGen/pecotmr")
```

If you have cloned the repository locally, you can install the package
with the `install_local` function from devtools. Assuming your working
directory contains the `pecotmr` repository, run this code to
install the package:

```R
devtools::install_local("pecotmr",build_vignettes = TRUE)
```

`devtools` should automatically install all required packages if
they are not installed already. **If automatic installation fails for some packages and you are unable to fix them, please try to use our conda package instead.**

To install using `conda` or `micromamba` for our released version,

```
micromamba install -c dnachun -c conda-forge -c bioconda r-pecotmr
```

To still use `devtools` to build `pecotmr` from source we suggest you use [pixi](https://pixi.sh) to set up a local environment with all dependencies:

```
.github/workflows/create_toml_from_yaml.sh .
pixi install
pixi shell
```

Now you should be able to build `pecotmr` from source.

## Developer's notes

+ When any changes are made to `roxygen2` markup or the C++ code in
the src directory, run `devtools::document()` to update the cpp11
registration code ([src/cpp11.cpp](src/cpp11.cpp) and
[R/cpp11.R](R/cpp11.R)), the package namespaces (see
[NAMESPACE](NAMESPACE)), and the package documentation files (in the
"man" subdirectory),

+ These are the R commands to build the website (make sure you are
connected to Internet while running these commands):

   ```R
   pkgdown::build_site_github_pages(new_process = FALSE, install = TRUE)
   ```

   The website is automatically deployed to https://statfungen.github.io/pecotmr
   on push to `main` via GitHub Actions.

+ To format R codes in the `R` folder,

   ```bash
   for i in `ls R/*.R`; do bash inst/misc/format_r_code.sh $i; done
   ```

+ After editing C++ code in the `src` directory, please use
[uncrustify][uncrustify] (version >=0.74.0, available from conda-forge) 
to format the code using configuration file
`inst/misc/uncrustify_default.cfg`. For example:

   ```bash
   uncrustify -c inst/misc/uncrustify_default.cfg --replace --no-backup -l CPP src/qtl_enrichment.cpp
   uncrustify -c inst/misc/uncrustify_default.cfg --replace --no-backup -l CPP src/qtl_enrichment.hpp
   ```

+ Prior to submitting the package to Bioconductor, the following modifications
need to be made: (1) remove the `Remotes:` entry in `DESCRIPTION`; and
(2) remove the `fastenloc.Rmd` vignette.

[BiocManager]: https://github.com/Bioconductor/BiocManager
[uncrustify]: https://github.com/uncrustify/uncrustify
