# `magicann`: Unveil the Magic of AnnData in R

`{magicann}` is an R package that brings the power and flexibility of AnnData to the R ecosystem, allowing you to effortlessly manipulate and analyze your single-cell data. With a touch of magic, `{magicann}` lets you work with backed h5ad and zarr files, directly access various slots (e.g. X, obs, var, obsm, obsp), or convert the data into SingleCellExperiment and Seurat objects. Our magician mascot ensures that your data analysis is always enchanting and fun!

## Features

* Seamless Integration: magicann bridges the gap between Python's AnnData and popular R single-cell analysis tools like SingleCellExperiment and Seurat.
* Efficient Data Handling: The package supports backed h5ad and zarr file formats, allowing you to work with large datasets efficiently.
* Flexible Data Access: Access obs, var, obsp, and other slots directly, giving you the freedom to manipulate your data as needed.
* Smooth Data Conversion: Effortlessly convert your AnnData objects to SingleCellExperiment or Seurat objects for downstream analysis.

## Design docs

* [Design document](doc/design.md): interface, OO-framework, approach
* [Challenges in reading h5ad files in R](doc/challenges.md)

## Installation

You can install the latest version of `{magicann}` from GitHub:

```r
# Install the devtools package if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install magicann from GitHub
devtools::install_github("scverse/magicann")
```

## Usage

Here's a quick example of how to use `{magicann}`:

```r
library(magicann)

# Read an h5ad file
adata <- read_h5ad("path/to/your/data.h5ad")

# Access the obs DataFrame
adata_obs <- adata$obs

# Convert the AnnData object to a SingleCellExperiment object
sce <- adata$to_sce()

# Convert the AnnData object to a Seurat object
seurat <- adata$to_seurat()
```
