#' Convert an `AnnData` to an `InMemoryAnnData`
#'
#' Convert another `AnnData` object to an [`InMemoryAnnData`] object
#'
#' @param adata An `AnnData` object to be converted to [`InMemoryAnnData`]
#'
#' @return An [`InMemoryAnnData`] object with the same data as the input `AnnData`
#'   object
#' @name as_InMemoryAnnData
#'
#' @family object converters
#'
#' @examples
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' ad$as_InMemoryAnnData()
NULL

#' Convert an `AnnData` to an `HDF5AnnData`
#'
#' Convert another `AnnData` object to an [`HDF5AnnData`] object
#'
#' @param adata An `AnnData` object to be converted to [`HDF5AnnData`]
#' @param file The file name (character) of the `.h5ad` file
#' @param compression The compression algorithm to use when writing the
#'  HDF5 file. Can be one of `"none"`, `"gzip"` or `"lzf"`. Defaults to
#' `"none"`.
#' @param mode The mode to open the HDF5 file:
#'
#'   * `a` creates a new file or opens an existing one for read/write
#'   * `r` opens an existing file for reading
#'   * `r+` opens an existing file for read/write
#'   * `w` creates a file, truncating any existing ones
#'   * `w-`/`x` are synonyms, creating a file and failing if it already exists
#'
#' @return An [`HDF5AnnData`] object with the same data as the input `AnnData`
#'   object.
#' @name as_HDF5AnnData
#'
#' @family object converters
#'
#' @examples
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5),
#' )
#' ad$as_HDF5AnnData("test.h5ad")
#' # remove file
#' file.remove("test.h5ad")
NULL

#' Convert an `AnnData` to a `SingleCellExperiment`
#'
#' Convert an `AnnData` object to a `SingleCellExperiment` object
#'
#' @param adata The `AnnData` object to convert
#' @param x_mapping A string specifying the name of the assay in the resulting
#'   `SingleCellExperiment` where the data in the `X` slot of `adata` will be
#'   mapped to
#' @param assays_mapping A named vector where names are names of `assays` in the
#'   resulting `SingleCellExperiment` object and values are keys of `layers` in
#'   `adata`. See below for details.
#' @param colData_mapping A named vector where names are columns of `colData` in
#'   the resulting `SingleCellExperiment` object and values are columns of `obs`
#'   in `adata`. See below for details.
#' @param rowData_mapping A named vector where names are columns of `rowData` in
#'   the resulting `SingleCellExperiment` object and values are columns of `var`
#'   in `adata`. See below for details.
#' @param reducedDims_mapping A named vector where names are names of
#'   `reducedDims` in the resulting `SingleCellExperiment` object and values are
#'   keys of `obsm` in `adata`. Alternatively, a named list where names are
#'   names of `reducedDims` in the resulting `SingleCellExperiment` object and
#'   values are vectors with the items `"sampleFactors"` and `"featureLoadings"`
#'   and/or `"metadata"`. See below for details.
#' @param colPairs_mapping A named vector where names are names of `colPairs` in
#'   the resulting `SingleCellExperiment` object and values are keys of `obsp`
#'   in `adata`. See below for details.
#' @param rowPairs_mapping A named vector where names are names of `rowPairs` in
#'   the resulting `SingleCellExperiment` object and values are keys of `varp`
#'   in `adata`. See below for details.
#' @param metadata_mapping A named vector where names are names of `metadata` in
#'   the resulting `SingleCellExperiment` object and values are keys of `uns` in
#'   `adata`. See below for details.
#'
#' @details
#'
#' ## Mapping arguments
#'
#' All mapping arguments expect a named character vector where names are the
#' names of the slot in the `SingleCellExperiment` object and values are the
#' keys of the corresponding slot of `adata`. If `NULL`, the conversion function
#' will guess which items to copy as described in the conversion tables
#' conversion table below. In most cases, the default is to copy all items using
#' the same names except where the correspondence between objects is unclear.
#' The `reducedDims_mapping` argument can also accept a more complex list format,
#' see below for details. To avoid copying anything to a slot, provide an empty
#' vector. If an unnamed vector is provided, the values will be used as names
#'
#' ### Examples:
#'
#'   - `NULL` will guess which items to copy as described in the conversion
#' table
#'   - `c(sce_item = "adata_item")` will copy `adata_item` from the slot in `adata` to
#' `sce_item` in the corresponding slot of new `SingleCellExperiment` object
#'   - `c()` will avoid copying anything to the slot
#'   - `c("adata_item")` is equivalent to `c(adata_item = "adata_item")`
#'
#' ## Conversion table
#'
# nolint start: line_length_linter
#'
#' | **From `AnnData`** | **To `SingleCellExperiment`** | **Example mapping argument** | **Default if `NULL`** |
#' |--------------------|-------------------------------|------------------------------|-----------------------|
#' | `adata$X` | `assays(sce)` | `x_mapping = "counts"` | The data in `adata$X` is copied to the assay named `X` |
#' | `adata$layers` | `assays(sce)` | `assays_mapping = c(counts = "counts")` | All items are copied by name |
#' | `adata$obs` | `colData(sce)` | `colData_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
#' | `adata$var` | `rowData(sce)` | `rowData_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
#' | `adata$obsm` | `reducedDims(sce)` | `reducedDims_mapping = c(pca = "X_pca")` **OR** `reducedDims_mapping = list(pca = c(obsm = "X_pca", varm = "PCs", uns = "pca_metadata"))`  | All items are copied by name without loadings except for `"X_pca"` for which loadings are added from `"PCs"` |
#' | `adata$obsp` | `colPairs(sce)` | `colPairs_mapping = c(nn = "connectivities")` | All items are copied by name |
#' | `adata$varp` | `rowPairs(sce)` | `rowPairs_mapping = c(gene_overlaps = "similarities")` | All items are copied by name |
#' | `adata$uns` | `metadata(sce)` | `uns_mapping = c(project_metadata = "metadata")` | All items are copied by name |
#'
# nolint end: line_length_linter
#'
#' ## The `reducedDims_mapping` argument
#'
#' For the simpler named vector format, the names should be the names of
#' `reducedDims` in the resulting `SingleCellExperiment` object, and the values
#' should be the keys of `obsm` in `adata`.
#'
#' For more advanced mapping, use the list format where each item is a vector
#' with the following names used to create a
#' [`SingleCellExperiment::LinearEmbeddingMatrix`] (if `featureLoadings` or
#' `metadata` is provided):
#'
#' - `sampleFactors`: a key of the `obsm` slot in `adata`,
#'   `adata$obsm[[sampleFactors]]` is passed to the `sampleFactors` argument
#' - `featureLoadings`: a key of the `varm` slot in `adata` (optional),
#'   `adata$varm[[featureLoadings]]` is passed to the `featureLoadings` argument
#' - `metadata`: a key of the `uns` slot in `adata` (optional),
#'   `adata$uns[[metadata]]` is passed to the `metadata` argument
#'
#' ## The `x_mapping` and `assays_mapping` arguments
#'
#' In order to specify where the data in `adata$X` will be stored in the
#' `assays(sce)` slot of the resulting object, you can use either the `x_mapping`
#' argument or the `assays_mapping` argument.
#' If you use `x_mapping`, it should be a string specifying the name of the layer
#' in `assays(sce)` where the data in `adata$X` will be stored.
#' If you use `assays_mapping`, it should be a named vector where names are names
#' of `assays(sce)` and values are keys of `layers` in `adata`.
#' In order to indicate the `adata$X` slot, you use `NA` as the value in the vector.
#' You may not specify both `NA` in `assays_mapping` and a value in `x_mapping`.
#' The name you provide for `x_mapping` may not be a name in `assays_mapping`.
#' You must provide an assay named `counts` or `data` in either `x_mapping` or
#' `assays_mapping`.
#'
#' @return A `SingleCellExperiment` object containing the requested data from
#'   `adata`
#' @name as_SingleCellExperiment
#'
#' @family object converters
#'
#' @examplesIf rlang::is_installed("SingleCellExperiment")
#'   ad <- AnnData(
#'     X = matrix(1:5, 3L, 5L),
#'     layers = list(A = matrix(5:1, 3L, 5L), B = matrix(letters[1:5], 3L, 5L)),
#'     obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'     var = data.frame(row.names = letters[1:5], gene = 1:5)
#'   )
#'
#'   # Default usage
#'   sce <- ad$as_SingleCellExperiment(
#'     assays_mapping = NULL,
#'     colData_mapping = NULL,
#'     rowData_mapping = NULL,
#'     reducedDims_mapping = NULL,
#'     colPairs_mapping = NULL,
#'     rowPairs_mapping = NULL,
#'     metadata_mapping = NULL
#'  )
NULL

#' Convert an `AnnData` to a `Seurat`
#'
#' Convert an `AnnData` object to a `Seurat` object
#'
#' @param adata The `AnnData` object to convert.
#' @param assay_name Name of the assay to be created in the new `Seurat` object
#' @param x_mapping A string specifying the name of the layer in the resulting
#'   `Seurat` object where the data in the `X` slot of `adata` will be mapped to
#' @param layers_mapping A named vector where names are names of `Layers` in the
#'   resulting `Seurat` object and values are keys of `layers` in `adata`. See
#'   below for details.
#' @param object_metadata_mapping A named vector where names are cell metadata
#'   columns in the resulting `Seurat` object and values are columns of `obs` in
#'   `adata`. See below for details.
#' @param assay_metadata_mapping A named vector where names are variable
#'   metadata columns in the assay of the resulting `Seurat` object and values
#'   are columns of `var` in `adata`. See below for details.
#' @param reduction_mapping A named vector where names are names of `Embeddings`
#'   in the resulting `Seurat` object and values are keys of `obsm` in `adata`.
#'   Alternatively, a named list where names are names of `Embeddings` in the
#'   resulting `Seurat` object and values are vectors with the items `"key"`,
#'   `"embeddings"` and (optionally) `"loadings"`. See below for details.
#' @param graph_mapping  A named vector where names are names of `Graphs` in the
#'   resulting `Seurat` object and values are keys of `obsp` in `adata`. See
#'   below for details.
#' @param misc_mapping A named vector where names are names of `Misc` in the
#'   resulting `Seurat` object and values are keys of `uns` in `adata`. See
#'   below for details.
#'
#' @details
#'
#' ## Mapping arguments
#'
#' All mapping arguments expect a named character vector where names are the
#' names of the slot in the `Seurat` object and values are the keys of the
#' corresponding slot of `adata`. If `NULL`, the conversion function will guess
#' which items to copy as described in the conversion tables conversion table
#' below. In most cases, the default is to copy all items using the same names
#' except where the correspondence between objects is unclear. The
#' `reduction_mapping` argument can also accept a more complex list format, see
#' below for details. To avoid copying anything to a slot, provide an empty
#' vector. If an unnamed vector is provided, the values will be used as names
#'
#' ### Examples:
#'
#' - `NULL` will guess which items to copy as described in the conversion table
#' - `c(seurat_item = "adata_item")` will copy `adata_item` from the slot in
#'   `adata` to `seurat_item` in the corresponding slot of new `Seurat` object
#' - `c()` will avoid copying anything to the slot
#' - `c("adata_item")` is equivalent to `c(adata_item = "adata_item")`
#'
#' ## Conversion table
#'
# nolint start: line_length_linter
#'
#' | **From `AnnData`** | **To `Seurat`** | **Example mapping argument** | **Default if `NULL`** |
#' |--------------------|-------------------------------|------------------------------|-----------------------|
#' | `adata$X` | `Layers(seurat)` | `x_mapping = "counts"` | The data in `adata$X` is copied to the assay named `X` |
#' | `adata$layers` | `Layers(seurat)` | `layers_mapping = c(counts = "counts")` | All items are copied by name |
#' | `adata$obs` | `seurat[[]]` | `object_metadata_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
#' | `adata$var` | `seurat[[assay_name]][[]]` | `assay_metadata_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
#' | `adata$obsm` | `Embeddings(sce)` | `reduction_mapping = c(pca = "X_pca")` **OR** `reduction_mapping = list(pca = c(key = "PC_", obsm = "X_pca", varm = "PCs"))`  | All items that can be coerced to a numeric matrix are copied by name without loadings except for `"X_pca"` for which loadings are added from `"PCs"` |
#' | `adata$obsp` | `Graphs(seurat)` | `graph_mapping = c(nn = "connectivities")` | All items are copied by name |
#' | `adata$varp` | _NA_  | _NA_ | There is no corresponding slot for `varp` |
#' | `adata$uns` | `Misc(seurat)` | `misc_mapping = c(project_metadata = "metadata")` | All items are copied by name |
#'
# nolint end: line_length_linter
#'
#' ## The `reduction_mapping` argument
#'
#' For the simpler named vector format, the names should be the names of
#' `Embeddings` in the resulting `Seurat` object, and the values
#' should be the keys of `obsm` in `adata`. A key will created from the `obsm`
#' key.
#'
#' For more advanced mapping, use the list format where each item is a vector
#' with the following names defining arguments to
#' [SeuratObject::CreateDimReducObject()]:
#'
#' - `key`: the key of the resulting [`SeuratObject::DimReduc`] object, passed
#'   to the `key` argument after processing with [SeuratObject::Key()]
#' - `embeddings`: a key of the `obsm` slot in `adata`,
#'   `adata$obsm[[embeddings]]` is passed to the `embeddings` argument
#' - `loadings`: a key of the `varm` slot in `adata` (optional),
#'   `adata$varm[[loadings]]` is passed to the `loadings` argument
#'
#' ## The `x_mapping` and `layers_mapping` arguments
#'
#' In order to specify where the data in `adata$X` will be stored in the
#' `Layers(seurat)` slot of the resulting object, you can use either the `x_mapping`
#' argument or the `layers_mapping` argument.
#' If you use `x_mapping`, it should be a string specifying the name of the layer
#' in `Layers(seurat)` where the data in `adata$X` will be stored.
#' If you use `layers_mapping`, it should be a named vector where names are names
#' of `Layers(seurat)` and values are keys of `layers` in `adata`.
#' In order to indicate the `adata$X` slot, you use `NA` as the value in the vector.
#' You may not specify both `NA` in `layers_mapping` and a value in `x_mapping`.
#' The name you provide for `x_mapping` may not be a name in `layers_mapping`.
#' You must provide a layer named `counts` or `data` in either `x_mapping` or
#' `layers_mapping`.
#'
#' @return A `Seurat` object containing the requested data from `adata`
#' @name as_Seurat
#'
#' @family object converters
#'
#' @examplesIf rlang::is_installed("Seurat")
#'   ad <- AnnData(
#'     X = matrix(1:5, 3L, 5L),
#'     obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'     var = data.frame(row.names = letters[1:5], gene = 1:5)
#'   )
#'
#'   # Default usage
#'   seurat <- ad$as_Seurat(
#'     assay_name = "RNA",
#'     x_mapping = "counts",
#'     layers_mapping = NULL,
#'     object_metadata_mapping = NULL,
#'     assay_metadata_mapping = NULL,
#'     reduction_mapping = NULL,
#'     graph_mapping = NULL,
#'     misc_mapping = NULL
#'   )
NULL
