#' Create an in-memory AnnData object.
#'
#' For more information on the functionality of an AnnData object, see [AnnData-usage]
#'
#' @param X See the `X` slot in [AnnData-usage]
#' @param layers See the `layers` slot in [AnnData-usage]
#' @param obs See the `obs` slot in [AnnData-usage]
#' @param var See the `var` slot in [AnnData-usage]
#' @param obsm See the `obsm` slot in [AnnData-usage]
#' @param varm See the `varm` slot in [AnnData-usage]
#' @param obsp See the `obsp` slot in [AnnData-usage]
#' @param varp See the `varp` slot in [AnnData-usage]
#' @param uns See the `uns` slot in [AnnData-usage]
#' @param shape Shape tuple (e.g. `c(n_obs, n_vars)`). Can be provided if both
#'   `X` or `obs` and `var` are not provided.
#'
#' @return An [InMemoryAnnData] object
#' @export
#'
#' @seealso [AnnData-usage] for details of `AnnData` structure and usage
#'
#' @examples
#' adata <- AnnData(
#'   X = matrix(1:12, nrow = 3, ncol = 4),
#'   obs = data.frame(
#'     row.names = paste0("obs", 1:3),
#'     n_counts = c(1, 2, 3),
#'     n_cells = c(1, 2, 3)
#'   ),
#'   var = data.frame(
#'     row.names = paste0("var", 1:4),
#'     n_cells = c(1, 2, 3, 4)
#'   )
#' )
#'
#' adata
AnnData <- function(
  X = NULL,
  obs = NULL,
  var = NULL,
  layers = NULL,
  obsm = NULL,
  varm = NULL,
  obsp = NULL,
  varp = NULL,
  uns = NULL,
  shape = NULL
) {
  InMemoryAnnData$new(
    X = X,
    obs = obs,
    var = var,
    layers = layers,
    obsm = obsm,
    varm = varm,
    obsp = obsp,
    varp = varp,
    uns = uns,
    shape = shape
  )
}

#' Convert to an `AnnData` object
#'
#' Convert other objects to an `AnnData` object. See the sections below for
#' details on how slots are mapped between objects. For more information on the
#' functionality of an `AnnData` object, see [AnnData-usage].
#'
#' @param x The object to convert
#' @param x_mapping A string specifying the data to map to the `X` slot. If
#'   `NULL`, no data will be copied to the `X` slot.
#' @param layers_mapping A named character vector where the names are keys of
#'   `layers` in the new `AnnData` object and values are the names of items in
#'   the corresponding slot of `x`. See below for default if `NULL` depending on
#'   the class of `x`.
#' @param obs_mapping A named character vector where the names are names of
#'   `obs` columns in the new `AnnData` object and values are the names of
#'   columns in the corresponding slot of `x`. See below for default if `NULL`
#'   depending on the class of `x`.
#' @param var_mapping A named character vector where the names are names of
#'   `var` columns in the new `AnnData` object and values are the names of
#'   columns in the corresponding slot of `x`. See below for default if `NULL`
#'   depending on the class of `x`.
#' @param obsm_mapping A named character vector where the names are keys of
#'   `obsm` in the new `AnnData` object and values are the names of items in the
#'   corresponding slot of `x`. See below for default if `NULL` depending on the
#'   class of `x`.
#' @param varm_mapping A named character vector where the names are keys of
#'   `varm` in the new `AnnData` object and values are the names of items in the
#'   corresponding slot of `x`. See below for default if `NULL` depending on the
#'   class of `x`.
#' @param obsp_mapping A named character vector where the names are keys of
#'   `obsp` in the new `AnnData` object and values are the names of items in the
#'   corresponding slot of `x`. See below for default if `NULL` depending on the
#'   class of `x`.
#' @param varp_mapping A named character vector where the names are keys of
#'   `varp` in the new `AnnData` object and values are the names of items in the
#'   corresponding slot of `x`. See below for default if `NULL` depending on the
#'   class of `x`.
#' @param uns_mapping A named character vector where the names are keys of `uns`
#'   in the new `AnnData` object and values are the names of items in the
#'   corresponding slot of `x`. See below for default if `NULL` depending on the
#'   class of `x`.
#' @param output_class The `AnnData` class to convert to. Must be one of
#'   `"HDF5AnnData"` or `"InMemoryAnnData"`.
#' @param ... Additional arguments passed to the generator function for
#'   `output_class`
#'
#' @section Details of mapping arguments:
#'
#'   All mapping arguments excepts for `x_mapping` expect a named character
#'   vector where names are the keys of the slot in the `AnnData` object and
#'   values are the names of items in the corresponding slot of `x`. If `NULL`,
#'   the conversion function will guess which items to copy as described in the
#'   conversion tables for each object type. In most cases, the default is to
#'   copy all items using the same names except where the correspondence between
#'   objects is unclear. To avoid copying anything to a slot, provide an empty
#'   vector. If an unnamed vector is provided, the values will be used as names
#'
#'   ## Examples:
#'
#'   - `NULL` will guess which items to copy as described in the conversion
#'     tables for each object type
#'   - `c(adata_item = "x_item")` will copy `x_item` from the slot in `x` to
#'     `adata_item` in the corresponding slot of new `AnnData` object
#'   - `c()` will avoid copying anything to the slot
#'   - `c("x_item")` is equivalent to `c(x_item = "x_item")`
#'
#' @section Converting from a `SingleCellExperiment` object:
#'
#'   This table describes how slots in a
#'   [`SingleCellExperiment::SingleCellExperiment`] object to the new `AnnData`
#'   object.
#'
#'   | **From `SingleCellExperiment`** | **To `AnnData`** | **Example mapping argument** | **Default if `NULL`** |
#'   |---------------------------------|------------------|------------------------------|-----------------------|
#'   | `assays(x)` | `adata$X` | `x_mapping = "counts"` | Nothing is copied to `X` |
#'   | `assays(x)` | `adata$layers` | `layers_mapping = c(counts = "counts")` | All items are copied by name |
#'   | `colData(x)` | `adata$obs` | `obs_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
#'   | `rowData(x)` | `adata$var` | `var_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
#'   | `reducedDims(x)` | `adata$obsm` | `obsm_mapping = c(X_pca = "pca")` | All items are copied by name |
#'   | `featureLoadings(reducedDims(x))` | `adata$varm` | `varm_mapping = c(PCs = "pca")` | Feature loadings from all [`SingleCellExperiment::LinearEmbeddingMatrix`] objects in `reducedDims(x)` |
#'   | `colPairs(x)` | `adata$obsp` | `obsp_mapping = c(connectivities = "RNA_nn")` | All items are copied by name |
#'   | `rowPairs(x)` | `adata$varp` | `varp_mapping = c(similarities = "gene_overlaps")` | All items are copied by name |
#'   | `metadata(x)` | `adata$uns` | `uns_mapping = c(metadata = "project_metadata")` | All items are copied by name |
#'
#' @section Converting from a `Seurat` object:
#'
#'   Only one assay can be converted from a [`SeuratObject::Seurat`] object to
#'   an `AnnData` object at a time. This can be controlled using the
#'   `assay_name` argument. By default, the current default assay will be used.
#'
#'   This table describes how slots in a [`SeuratObject::Seurat`] object to the
#'   new `AnnData` object.
#'
#'   | **From `Seurat`** | **To `AnnData`** | **Example mapping argument** | **Default if `NULL`** |
#'   |-------------------|------------------|------------------------------|-----------------------|
#'   | `Layers(x)` | `adata$X` | `x_mapping = "counts"` | Nothing is copied to `X` |
#'   | `Layers(x)` | `adata$layers` | `layers_mapping = c(counts = "counts")` | All items are copied by name |
#'   | `Misc(x)` | `adata$obs` | `obs_mapping = c(n_counts = "n_counts", cell_type = "CellType")` | All columns are copied by name |
#'   | `Misc(x[[assay_name]])` | `adata$var` | `var_mapping = c(n_cells = "n_cells", pct_zero = "PctZero")` | All columns are copied by name |
#'   | `Embeddings(x)` | `adata$obsm` | `obsm_mapping = c(X_pca = "pca")` | All embeddings matching `assay_name` are copied by name |
#'   | `Loadings(x)` | `adata$varm` | `varm_mapping = c(PCs = "pca")` | All valid loadings are copied by name |
#'   | `Graphs(x)` | `adata$obsp` | `obsp_mapping = c(connectivities = "RNA_nn")` | All graphs matching `assay_name` are copied by name |
#'   | `Misc(x)` | `adata$varp` | `varp_mapping = c(similarities = "gene_overlaps")` | No data is copied to `varp` |
#'   | `Misc(x)` | `adata$uns` | `uns_mapping = c(metadata = "project_metadata")` | All items are copied by name |
#'
#' @return An `AnnData` object of the class requested by `output_class`
#'   containing the data specified in the mapping arguments.
#'
#' @export
#'
#' @examplesIf rlang::is_installed("Seurat")
#' # Convert a Seurat object to an AnnData object
#' library(Seurat)
#'
#' counts <- matrix(rbinom(20000, 1000, .001), nrow = 100)
#' obj <- CreateSeuratObject(counts = counts)
#' obj <- NormalizeData(obj)
#' obj <- FindVariableFeatures(obj)
#' obj <- ScaleData(obj)
#' obj <- RunPCA(obj, npcs = 10L)
#' obj <- FindNeighbors(obj)
#' obj <- RunUMAP(obj, dims = 1:10)
#'
#' as_AnnData(obj)
#'
#' @examplesIf rlang::is_installed("SingleCellExperiment")
#' # Convert a SingleCellExperiment object to an AnnData object
#' library(SingleCellExperiment)
#'
#' sce <- SingleCellExperiment(
#'   assays = list(counts = matrix(1:5, 5L, 3L)),
#'   colData = DataFrame(cell = 1:3, row.names = paste0("Cell", 1:3)),
#'   rowData = DataFrame(gene = 1:5, row.names = paste0("Gene", 1:5))
#' )
#'
#' as_AnnData(sce)
as_AnnData <- function(
  x,
  x_mapping = NULL,
  layers_mapping = NULL,
  obs_mapping = NULL,
  var_mapping = NULL,
  obsm_mapping = NULL,
  varm_mapping = NULL,
  obsp_mapping = NULL,
  varp_mapping = NULL,
  uns_mapping = NULL,
  output_class = c("InMemory", "HDF5AnnData"),
  ...
) {
  UseMethod("as_AnnData", x)
}

#' @rdname as_AnnData
#' @export
as_AnnData.SingleCellExperiment <- function(
    x,
    x_mapping = NULL,
    layers_mapping = NULL,
    obs_mapping = NULL,
    var_mapping = NULL,
    obsm_mapping = NULL,
    varm_mapping = NULL,
    obsp_mapping = NULL,
    varp_mapping = NULL,
    uns_mapping = NULL,
    output_class = c("InMemory", "HDF5AnnData"),
    ...
) {
  from_SingleCellExperiment(
    sce = x,
    x_mapping = x_mapping,
    layers_mapping = layers_mapping,
    obs_mapping = obs_mapping,
    var_mapping = var_mapping,
    obsm_mapping = obsm_mapping,
    varm_mapping = varm_mapping,
    obsp_mapping = obsp_mapping,
    varp_mapping = varp_mapping,
    uns_mapping = uns_mapping,
    output_class = output_class,
    ...
  )
}

#' @param assay_name For [`SeuratObject::Seurat`] objects, the name of the assay
#'   to be converted. If `NULL`, the default assay will be used
#'   ([SeuratObject::DefaultAssay()]). This is ignored for other objects.
#'
#' @rdname as_AnnData
#' @export
as_AnnData.Seurat <- function(
    x,
    assay_name = NULL,
    x_mapping = NULL,
    layers_mapping = NULL,
    obs_mapping = NULL,
    var_mapping = NULL,
    obsm_mapping = NULL,
    varm_mapping = NULL,
    obsp_mapping = NULL,
    varp_mapping = NULL,
    uns_mapping = NULL,
    output_class = c("InMemory", "HDF5AnnData"),
    ...
) {
  from_Seurat(
    seurat_object = x,
    assay_name = assay_name,
    x_mapping = x_mapping,
    layers_mapping = layers_mapping,
    obs_mapping = obs_mapping,
    var_mapping = var_mapping,
    obsm_mapping = obsm_mapping,
    varm_mapping = varm_mapping,
    obsp_mapping = obsp_mapping,
    varp_mapping = varp_mapping,
    uns_mapping = uns_mapping,
    output_class = output_class,
    ...
  )
}
