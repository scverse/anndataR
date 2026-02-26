# Write H5AD

Write an H5AD file

## Usage

``` r
write_h5ad(
  object,
  path,
  compression = c("none", "gzip", "lzf"),
  chunk_size = "auto",
  mode = c("w-", "r", "r+", "a", "w", "x"),
  ...
)
```

## Arguments

- object:

  The object to write, either a
  [`SingleCellExperiment::SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  or a
  [`SeuratObject::Seurat`](https://satijalab.github.io/seurat-object/reference/Seurat-class.html)
  object

- path:

  Path of the file to write to

- compression:

  The compression algorithm to use when writing the HDF5 file. Can be
  one of `"none"`, `"gzip"` or `"lzf"`. Defaults to `"none"`.

- chunk_size:

  The target chunk size in bytes to use when writing HDF5 datasets. When
  `"auto"` (default), the chunk size is determined automatically using
  an algorithm that mimics h5py's auto-chunking behaviour. Set to `NULL`
  to disable chunking (contiguous storage, the rhdf5 default), or a
  number to use a specific target size in bytes. This only affects
  array-like datasets; scalar values are unaffected.

- mode:

  The mode to open the HDF5 file.

  - `a` creates a new file or opens an existing one for read/write

  - `r+` opens an existing file for read/write

  - `w` creates a file, truncating any existing ones

  - `w-`/`x` are synonyms creating a file and failing if it already
    exists

- ...:

  Additional arguments passed to
  [`as_AnnData()`](https://anndataR.scverse.org/reference/as_AnnData.md)

## Value

`path` invisibly

## Details

### Compression

Compression is currently not supported for Boolean arrays, they will be
written uncompressed.

### `NULL` values

For compatibility with changes in Python **anndata** 0.12.0, `NULL`
values in `uns` are written to H5AD files as a `NULL` dataset (instead
of not being written at all). To disable this behaviour, set
`option(anndataR.write_null = FALSE)`. This may be required to allow the
file to be read by older versions of Python **anndata**.

## Examples

``` r
adata <- AnnData(
  X = matrix(1:5, 3L, 5L),
  layers = list(
    A = matrix(5:1, 3L, 5L),
    B = matrix(letters[1:5], 3L, 5L)
  ),
  obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
  var = data.frame(row.names = letters[1:5], gene = 1:5)
)
h5ad_file <- tempfile(fileext = ".h5ad")
adata$write_h5ad(h5ad_file)

# Write a SingleCellExperiment as an H5AD
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  ncells <- 100
  counts <- matrix(rpois(20000, 5), ncol = ncells)
  logcounts <- log2(counts + 1)

  pca <- matrix(runif(ncells * 5), ncells)
  tsne <- matrix(rnorm(ncells * 2), ncells)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts, logcounts = logcounts),
    reducedDims = list(PCA = pca, tSNE = tsne)
  )

  adata <- as_AnnData(sce)
  h5ad_file <- tempfile(fileext = ".h5ad")
  adata$write_h5ad(h5ad_file)
}

# Write a Seurat as a H5AD
if (requireNamespace("Seurat", quietly = TRUE)) {
  library(Seurat)

  counts <- matrix(1:15, 5L, 3L)
  dimnames(counts) <- list(
    LETTERS[1:5],
    letters[1:3]
  )
  cell.metadata <- data.frame(
    row.names = letters[1:3],
    cell = 1:3
  )
  obj <- CreateSeuratObject(counts, meta.data = cell.metadata)
  gene.metadata <- data.frame(
    row.names = LETTERS[1:5],
    gene = 1:5
  )
  obj[["RNA"]] <- AddMetaData(GetAssay(obj), gene.metadata)

  adata <- as_AnnData(obj)
  h5ad_file <- tempfile(fileext = ".h5ad")
  adata$write_h5ad(h5ad_file)
}
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
```
