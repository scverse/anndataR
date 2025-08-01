# {anndataR}: An R package for working with AnnData objects <img src="man/figures/logo.png" align="right" alt="anndataR logo" width=120 />
<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/scverse/anndataR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scverse/anndataR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**{anndataR}** aims to make the `AnnData` format a first-class citizen in the R ecosystem, and to make it easy to work with AnnData files in R, either directly or by converting them to a `SingleCellExperiment` or `Seurat` object.

**{anndataR}** is an scverse® community project maintained by [Data Intuitive](https://data-intuitive.com/), and is fiscally sponsored by the [Chan Zuckerberg Initiative](https://chanzuckerberg.com/).

## Features of {anndataR}

- Provide an `R6` class to work with `AnnData` objects in R (either in-memory or on-disk)
- Read/write `*.h5ad` files natively
- Convert to/from `SingleCellExperiment` objects
- Convert to/from `Seurat` objects

You can find the status of development of **{anndataR}** on the [feature tracking page](https://anndatar.data-intuitive.com/articles/design.html#feature-tracking) of the package website. 
Please [report](https://github.com/scverse/anndataR/issues) any issues you encounter.

## Installation

You can install **{anndataR}** from Bioconductor using **BiocManager**:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("anndataR")
```

Or you can install the development version of **{anndataR}** from GitHub like so:

``` r
# install.packages("pak")
pak::pak("scverse/anndataR")
```

You will need to install additional dependencies, depending on
the task you want to perform.

- To read/write `*.h5ad` files, install [rhdf5](https://www.bioconductor.org/packages/rhdf5):  
  `BiocManager::install("rhdf5")`
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

In order to browse these vignettes locally, you need to build them during installation:

```
options(pkg.build_vignettes = TRUE)
pak::pak("scverse/anndataR")
```

Take note that you need all suggested dependencies available, and that building them can take some time.

- [**Getting started**](https://anndatar.data-intuitive.com/articles/anndataR.html): An introduction to the package and its features.  
  `vignette("anndataR", package = "anndataR")`
- [**Read/write `Seurat` objects**](https://anndatar.data-intuitive.com/articles/usage_seurat.html): How to convert between `AnnData` and `Seurat` objects.  
  `vignette("usage_seurat", package = "anndataR")`
- [**Read/write `SingleCellExperiment` objects**](https://anndatar.data-intuitive.com/articles/usage_singlecellexperiment.html): How to convert between `AnnData` and `SingleCellExperiment` objects 
  `vignette("usage_singlecellexperiment", package = "anndataR")`
- [**Software Design**](https://anndatar.data-intuitive.com/articles/software_design.html): An overview of the design of the package
- [**Development Status**](https://anndatar.data-intuitive.com/articles/development_status.html): An overview of the development status of the package
- [**Known Issues**](https://anndatar.data-intuitive.com/articles/known_issues.html): An overview of known issues with the package.  
  `vignette("known_issues", package = "anndataR")`
