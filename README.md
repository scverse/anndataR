# {anndataR}: An R package for working with AnnData objects

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN status](https://www.r-pkg.org/badges/version/anndataR.png)](https://CRAN.R-project.org/package=anndataR)
[![R-CMD-check](https://github.com/scverse/anndataR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scverse/anndataR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**{anndataR}** aims to make the AnnData format a first-class citizen in
the R ecosystem, and to make it easy to work with AnnData files in R,
either directly or by converting them to a SingleCellExperiment or Seurat
object.

**{anndataR}** is an scverse community project maintained by [Data Intuitive](https://data-intuitive.com/), and is fiscally sponsored by the [Chan Zuckerberg Initiative](https://chanzuckerberg.com/).


## Features of {anndataR}

- Provide an `R6` class to work with AnnData objects in R (either in-memory or on-disk).
- Read/write `*.h5ad` files natively
- Convert to/from `SingleCellExperiment` objects
- Convert to/from `Seurat` objects

> [!WARNING]
>
> This package is still in the experimental stage, and may not work as
> expected. You can find the status of development of anndataR on the
> [feature tracking page](https://anndatar.data-intuitive.com/articles/design.html#feature-tracking)
> of the website. Please [report](https://github.com/scverse/anndataR/issues) any issues you encounter.

## Installation

You can install the development version of **{anndataR}** like so:

``` r
# install.packages("pak")
pak::pak("scverse/anndataR")
```

You will need to install additional dependencies, depending on
the task you want to perform.

- To read/write `*.h5ad` files, install [hdf5r](https://cran.r-project.org/package=hdf5r):  
  `install.packages("hdf5r")`
- To convert to/from `SingleCellExperiment` objects, install [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html):  
  `BiocManager::install("SingleCellExperiment")`
- To convert to/from `Seurat` objects, install [SeuratObject](https://cran.r-project.org/package=SeuratObject):  
  `install.packages("SeuratObject")`

Alternatively, you can install all suggested dependencies at once:

``` r
pak::pak("scverse/anndataR", dependencies = TRUE)
```

## Getting started

The best way to get started with **{anndataR}** is to explore the package vignettes (available at https://anndatar.data-intuitive.com/articles/).

- **Getting started**: An introduction to the package and its features.  
  `vignette("anndataR", package = "anndataR")`
- **Reading and writing H5AD files**: How to read and write `*.h5ad` files.  
  `vignette("usage_h5ad", package = "anndataR")`
- **Converting to/from Seurat objects**: How to convert between `AnnData` and `Seurat` objects.  
  `vignette("usage_seurat", package = "anndataR")`
- **Converting to/from SingleCellExperiment objects**: How to convert between `AnnData` and `SingleCellExperiment` objects.  
  `vignette("usage_singlecellexperiment", package = "anndataR")`
- **Software Design**: An overview of the design of the package.  
  `vignette("software_design", package = "anndataR")`
- **Development Status**: An overview of the development status of the package.  
  `vignette("development_status", package = "anndataR")`
- **Known Isses**: An overview of known issues with the package.  
  `vignette("known_issues", package = "anndataR")`
