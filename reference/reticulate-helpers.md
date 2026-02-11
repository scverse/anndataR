# Reticulate Helper Functions for AnnData Conversion

This file contains helper functions that enable seamless conversion
between R and Python AnnData objects using the reticulate package. These
functions provide automatic S3 method dispatch for converting AnnData
objects across the R-Python boundary.

## Usage

``` r
# S3 method for class 'collections.abc.Mapping'
py_to_r(x)

# S3 method for class 'anndata._core.anndata.AnnData'
py_to_r(x)

# S3 method for class 'AbstractAnnData'
r_to_py(x, convert = TRUE)
```

## Arguments

- x:

  An AbstractAnnData object (any anndataR implementation)

- convert:

  Whether to convert the result (passed to reticulate)

## Value

A
[`ReticulateAnnData`](https://anndataR.data-intuitive.com/reference/ReticulateAnnData.md)
object wrapping the Python object

A Python AnnData object

## Details

The main conversion functions include:

- `py_to_r.anndata._core.anndata.AnnData`: Converts Python AnnData
  objects to R
  [ReticulateAnnData](https://anndataR.data-intuitive.com/reference/ReticulateAnnData.md)
  objects

- `r_to_py.AbstractAnnData`: Converts R
  [AbstractAnnData](https://anndataR.data-intuitive.com/reference/AbstractAnnData.md)
  objects to Python AnnData objects

- `py_to_r.collections.abc.Mapping`: Converts Python mapping objects to
  R lists

These functions are automatically registered as S3 methods and are
called when using
[`reticulate::py_to_r()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
and
[`reticulate::r_to_py()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
on compatible objects.

## See also

Other object converters:
[`as_AnnData()`](https://anndataR.data-intuitive.com/reference/as_AnnData.md),
[`as_HDF5AnnData()`](https://anndataR.data-intuitive.com/reference/as_HDF5AnnData.md),
[`as_InMemoryAnnData()`](https://anndataR.data-intuitive.com/reference/as_InMemoryAnnData.md),
[`as_ReticulateAnnData()`](https://anndataR.data-intuitive.com/reference/as_ReticulateAnnData.md),
[`as_Seurat()`](https://anndataR.data-intuitive.com/reference/as_Seurat.md),
[`as_SingleCellExperiment()`](https://anndataR.data-intuitive.com/reference/as_SingleCellExperiment.md)

## Examples

``` r
# \donttest{
# Requires Python anndata to be installed
if (requireNamespace("reticulate", quietly = TRUE) &&
      reticulate::py_module_available("anndata")) {

  library(reticulate)

  # Create Python AnnData object
  ad_py <- import("anndata", convert = FALSE)
  py_adata <- ad_py$AnnData(X = r_to_py(matrix(1:12, 3, 4)))

  # Automatic conversion to R (uses py_to_r.anndata._core.anndata.AnnData)
  r_adata <- py_to_r(py_adata)

  # Automatic conversion back to Python (uses r_to_py.AbstractAnnData)
  py_adata2 <- r_to_py(r_adata)
}
# }
```
