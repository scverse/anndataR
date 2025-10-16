# Agents Instructions for anndataR

## Project Overview

**anndataR** is a Bioconductor R package that provides native R support for reading/writing `.h5ad` (AnnData) files and bidirectional conversion between AnnData objects and popular Bioconductor/Seurat formats.

**Key Features:**

- Native R reading/writing of `.h5ad` files (HDF5 backend)
- Three AnnData implementations: `InMemoryAnnData`, `HDF5AnnData`, `ReticulateAnnData`
- Bidirectional conversion: SingleCellExperiment ↔ AnnData ↔ Seurat
- S4 coercion system for seamless type conversion
- Python interoperability via reticulate

## Architecture

### R6 Class Hierarchy

    AbstractAnnData (abstract base class)
    ├── InMemoryAnnData (all data in RAM)
    ├── HDF5AnnData (HDF5-backed, lazy loading)
    ├── ReticulateAnnData (Python anndata wrapper via reticulate)
    └── AnnDataView (lazy subsetting view, no data copying)

**When to use each:**

- **InMemoryAnnData**: Default for most operations, fast for small/medium datasets
- **HDF5AnnData**: Large datasets exceeding RAM, on-disk persistence
- **ReticulateAnnData**: Interoperability with Python anndata library, testing against Python behavior
- **AnnDataView**: Lazy subsetting/slicing without copying data; convert to concrete implementation when needed

### Core Slots (AnnData spec)

- `X`: Main matrix (n_obs × n_vars)
- `layers`: Named list of alternative matrices with same dimensions as X
- `obs`: DataFrame with observation (cell) metadata (n_obs rows)
- `var`: DataFrame with variable (gene) metadata (n_vars rows)
- `obsm`: List of observation-aligned matrices (e.g., PCA, UMAP embeddings)
- `varm`: List of variable-aligned matrices
- `obsp`: List of observation pairwise matrices (e.g., cell-cell distances)
- `varp`: List of variable pairwise matrices
- `uns`: Unstructured metadata (arbitrary nested lists/dicts)

**Important validation rules:**

- All matrices in X/layers must have same shape: `[n_obs, n_vars]`
- Row names of obs/var define observation/variable names
- obsm/varm matrices must align: first dimension matches n_obs/n_vars
- obsp/varp must be square pairwise matrices

## Critical Implementation Details

### 1. obs_names/var_names Architecture

obs_names and var_names are now stored separately from obs/var data.frames:

```r
# InMemoryAnnData and HDF5AnnData store names separately
private$.obs_names  # Character vector of observation names
private$.var_names  # Character vector of variable names

# obs/var data.frames are stored WITHOUT rownames internally
private$.obs  # Data.frame with NULL rownames
private$.var  # Data.frame with NULL rownames

# Dimnames are added ON-THE-FLY when users access data
ad$obs        # Returns data.frame with rownames = obs_names
ad$X          # Returns matrix with dimnames from obs_names/var_names
```

**Why this matters:**

- All matrix data (X, layers, obsm, varm, obsp, varp) is stored internally **without dimnames**
- Dimnames are dynamically added via `.add_matrix_dimnames()` and `.add_obsvar_dimnames()` helper methods
- This ensures consistency with Python anndata where obs/var names are separate from the data
- When setting obs/var, extract rownames first if present, then strip them

**Pattern for setters:**

```r
# Extract names before validation
if (!is.null(value) && has_row_names(value)) {
  private$.obs_names <- rownames(value)
  rownames(value) <- NULL
}
private$.obs <- private$.validate_obsvar_dataframe(value, "obs")
```

**Use `obs_names`/`var_names` properties directly:**

- Don't use `rownames(ad$obs)` or `colnames(ad$X)` in internal code
- Use `ad$obs_names` and `ad$var_names` for reliable access to names

### 2. AnnDataView for Lazy Subsetting

`AnnDataView` provides lazy subsetting without data copying:

```r
# Create a view with S3 [ operator
ad <- AnnData(
  X = matrix(1:15, 3L, 5L),
  obs = data.frame(row.names = LETTERS[1:3], cell_type = c("A", "B", "A"))
)
view <- ad[ad$obs$cell_type == "A", ]  # Returns AnnDataView, no data copied

# Convert to concrete implementation when needed
result <- view$as_InMemoryAnnData()  # Now subsetting is applied
```

**Key characteristics:**

- Inherits from AbstractAnnData
- Stores base AnnData object and subset indices
- All getters apply subsetting on-the-fly
- Setters are **disabled** - must convert to concrete implementation first
- Use `.apply_subset()` helper for matrix/data.frame subsetting
- Use `.apply_vector_subset()` for vector subsetting (obs_names, var_names)

**Testing pattern:** See `tests/testthat/test-AnnDataView.R` for usage examples

### 3. S3 Methods for AbstractAnnData

Standard R S3 methods now work on all AnnData objects:

```r
# Dimension methods
dim(ad)           # [n_obs, n_vars]
nrow(ad)          # n_obs
ncol(ad)          # n_vars
dimnames(ad)      # list(obs_names, var_names)
rownames(ad)      # obs_names
colnames(ad)      # var_names

# Subsetting with [ operator
ad[1:5, ]         # Subset observations, returns AnnDataView
ad[, 1:10]        # Subset variables, returns AnnDataView
ad[1:5, 1:10]     # Subset both, returns AnnDataView
```

**Implementation:** See `R/AbstractAnnData-s3methods.R` for all S3 methods

### 4. HDF5 File Management

**Pattern:** HDF5AnnData manages file handles with lifecycle:

```r
# Automatic closure on finalization
adata <- read_h5ad("file.h5ad")  # close_on_finalize = TRUE

# Manual closure required for explicit construction
adata <- HDF5AnnData$new("file.h5ad")
adata$close()  # Must call explicitly

# Check validity before operations
private$.check_file_valid()  # Throws if handle closed
```

**Testing pattern:** `test-h5ad-fileclosure.R` validates handles close properly.

### 5. Validation Patterns

**Alignment validation:**

```r
# All setters use validation helpers from AbstractAnnData
private$.validate_aligned_array(
  value,
  "X",
  shape = c(self$n_obs(), self$n_vars()),
  expected_rownames = rownames(self),
  expected_colnames = colnames(self)
)
```

**Dimension checking:**

```r
private$.validate_aligned_mapping(
  value,
  "layers",
  c(self$n_obs(), self$n_vars()),
  expected_rownames = rownames(self),
  expected_colnames = colnames(self)
)
```

## Testing Infrastructure

### Conditional Test Skipping

**Pattern:** Helper functions check for optional dependencies:

```r
# tests/testthat/helper-skip_if_no_anndata.R
skip_if_no_anndata <- function() {
  if (!rlang::is_installed("reticulate")) {
    skip("reticulate not installed")
  }
  
  # Check Python anndata module
  if (!reticulate::py_module_available("anndata")) {
    skip("Python anndata not available")
  }
}

# Usage in tests
test_that("ReticulateAnnData works", {
  skip_if_no_anndata()
  # ... test code ...
})
```

**Available skip helpers:**

- `skip_if_no_anndata()`: Python anndata + reticulate
- `skip_if_no_dummy_anndata()`: anndataR.testdata Python module (for roundtrip tests)
- `skip_if_no_h5diff()`: h5diff CLI tool

### Roundtrip Testing Pattern

**Goal:** Ensure R read/write matches Python anndata behavior

```r
# tests/testthat/test-roundtrip-X.R example
test_that("Dense X roundtrips correctly", {
  skip_if_no_dummy_anndata()
  
  # Generate test data with known structure
  adata <- generate_dataset(
    X_type = "dense",
    obs = data.frame(row.names = LETTERS[1:10]),
    var = data.frame(row.names = letters[1:5])
  )
  
  # Write R → H5AD
  tmp <- tempfile(fileext = ".h5ad")
  adata$write_h5ad(tmp)
  
  # Read back and verify
  adata2 <- read_h5ad(tmp)
  expect_equal(adata$X, adata2$X)
  
  # Compare with Python via dummy_anndata
  py <- reticulate::import("anndataR.testdata")
  py_adata <- py$dummy_anndata(X = "dense")
  expect_equal_py(adata, py_adata)  # Custom helper
})
```

**Custom helpers:**

- `expect_equal_py(r_adata, py_adata)`: Compare R and Python AnnData objects
- `generate_dataset()`: Create test AnnData with configurable matrix types

### Mock Data Generation

**Usage:**

```r
# R/generate_dataset.R
adata <- generate_dataset(
  X_type = "dense",           # "dense", "sparse", "csparse", "rsparse"
  obs = data.frame(row.names = LETTERS[1:100]),
  var = data.frame(row.names = letters[1:50]),
  n_layers = 2,
  obs_has_row_names = TRUE    # Toggle for testing edge cases
)
```

## Conversion System

### From External Formats

**Pattern:** `from_*()` functions for explicit conversion:

```r
# from_SingleCellExperiment.R
adata <- from_SingleCellExperiment(
  sce,
  output_class = "InMemoryAnnData",  # or "HDF5AnnData"
  X_name = "counts",                  # Which assay → X
  layers = c("logcounts"),            # Additional assays → layers
  uns_keys = c("pca")                 # Which metadata → uns
)

# from_Seurat.R
adata <- from_Seurat(
  seurat_obj,
  output_class = "InMemoryAnnData",
  assay = "RNA"                       # Which assay to extract
)
```

**Guessing pattern:** `from_Seurat()` uses helper functions to intelligently map:

- `.from_Seurat_guess_layers()`: Identify assay data slots
- `.from_Seurat_guess_obsms()`: Extract dimensionality reductions (PCA, UMAP)
- `.from_Seurat_guess_uns()`: Map miscellaneous metadata

### To External Formats

**Pattern:** `as_*()` functions with optional parameters:

```r
# as_SingleCellExperiment.R
sce <- as_SingleCellExperiment(
  adata,
  X_name = "counts",        # What to call X assay
  layer_names = NULL        # Which layers to include (NULL = all)
)

# as_Seurat.R
seurat <- as_Seurat(
  adata,
  assay_name = "RNA"        # Assay name in Seurat object
)
```

## Development Workflow

### Build and Test Commands

```bash
# Full R CMD check (Bioconductor standards)
R CMD build .
R CMD check --as-cran anndataR_*.tar.gz

# Quick test suite
R -e 'devtools::test()'

# Run specific test file
R -e 'devtools::test(filter = "roundtrip-X")'

# Lint checking
R -e 'lintr::lint_package()'

# Check code formatting
air format --check .

# Reformat code
air format .
```

### Lintr Configuration

**Key rules:**

- Line length: 120 characters max
- Use `paste()` to wrap long strings in test helpers
- Prefer explicit `::` for package functions in examples

**Example fix:**

```r
# ❌ Too long
expect_warning(as(adata, "Seurat"), "Consider using as_Seurat() for more control over the conversion")

# ✅ Wrapped
expect_warning(
  as(adata, "Seurat"),
  paste(
    "Consider using as_Seurat() for more control",
    "over the conversion"
  )
)
```

### Dependency Management

**Pattern:** Check for optional packages before use:

```r
# R/check_requires.R
check_requires <- function(package, reason = NULL) {
  if (!rlang::is_installed(package)) {
    msg <- c(
      "!" = "Package {.pkg {package}} is required",
      "i" = "Install with: {.code install.packages(\"{package}\")}"
    )
    if (!is.null(reason)) {
      msg <- c(msg, "i" = reason)
    }
    cli::cli_abort(msg)
  }
}

# Usage
as_Seurat <- function(adata, ...) {
  check_requires("SeuratObject", "for converting to Seurat objects")
  # ... conversion code ...
}
```


## Common Pitfalls

### 1. obs_names/var_names Separation

**Problem:** Trying to access obs/var names via rownames instead of dedicated properties

**Wrong:**
```r
# ❌ Internal code should not use this pattern
obs_ids <- rownames(adata$obs)
gene_ids <- colnames(adata$X)
```

**Correct:**
```r
# ✅ Always use dedicated properties
obs_ids <- adata$obs_names
gene_ids <- adata$var_names
```

**Why:** obs_names/var_names are stored separately and added on-the-fly to user-facing data

### 2. Row Names Handling

**Problem:** R data.frames require unique row names, but AnnData allows duplicates

**Solution:** `generate_dataframe.R` has `obs_has_row_names` parameter for testing edge cases

### 3. Sparse Matrix Compatibility

**Problem:** Different sparse matrix classes (Matrix::dgCMatrix vs Matrix::dgRMatrix)

**Solution:** Use `generate_dataset(X_type = "csparse")` for CSC, `"rsparse"` for CSR testing

### 4. NULL vs Empty List in uns

**Problem:** Python distinguishes `None` from `{}`; R treats `NULL` differently

**Solution:** As of anndata 0.12.0, write `NULL` as empty HDF5 dataset. Controlled by `options(anndataR.write_null = TRUE)`

### 5. Closure Variable Capture in Factories

**Problem:** Loop variables not captured in closure

```r
# ❌ WRONG - all handlers reference final iteration value
for (class in classes) {
  setAs(class, "Seurat", function(from) convert(from, class))
}

# ✅ CORRECT - force() captures variable value
.make_convert_handler <- function(convert_fn, from_str, to_str) {
  force(convert_fn)
  force(from_str)
  force(to_str)
  function(from) convert_fn(from)
}
```

## File Organization

### Source Files by Category

**Core Classes:**

- `R/AbstractAnnData.R`: Base class with abstract slots and validation
- `R/AbstractAnnData-s3methods.R`: S3 methods (dim, nrow, ncol, dimnames, `[`)
- `R/InMemoryAnnData.R`: RAM-based implementation
- `R/HDF5AnnData.R`: HDF5-backed implementation
- `R/ReticulateAnnData.R`: Python wrapper (experimental)
- `R/AnnDataView.R`: Lazy view for subsetting without data copying

**I/O:**

- `R/read_h5ad.R`, `R/read_h5ad_helpers.R`: Reading HDF5 files
- `R/write_h5ad.R`, `R/write_h5ad_helpers.R`: Writing HDF5 files
- `R/write_hdf5_helpers.R`: Low-level HDF5 utilities

**Conversion:**

- `R/as-coercions.R`: S4 coercion registration
- `R/as_AnnData.R`: Generic converter to AnnData
- `R/as_SingleCellExperiment.R`: AnnData → SCE
- `R/as_Seurat.R`: AnnData → Seurat
- `R/from_SingleCellExperiment.R`: SCE → AnnData
- `R/from_Seurat.R`: Seurat → AnnData

**Testing:**

- `R/generate_dataset.R`, `R/generate_*.R`: Mock data generation
- `tests/testthat/helper-*.R`: Test utilities and skip helpers
- `tests/testthat/test-roundtrip-*.R`: Python compatibility tests
- `tests/testthat/test-as-*.R`: Conversion validation
- `tests/testthat/test-AnnDataView.R`: Lazy subsetting tests
- `tests/testthat/test-AbstractAnnData-s3methods.R`: S3 method tests

**Utilities:**

- `R/check_requires.R`: Dependency validation
- `R/ui.R`: CLI messaging helpers
- `R/utils.R`: Miscellaneous helpers
- `R/known_issues.R`: Track known bugs/limitations

## Documentation Standards

### Roxygen2 Patterns

**Class documentation:**

```r
#' @title InMemoryAnnData
#' @description Implementation of an in-memory AnnData object.
#' @seealso [AnnData-usage] for details on creating and using AnnData objects
#' @family AnnData classes
#' @examples
#' adata <- AnnData(X = matrix(1:15, 3L, 5L), ...)
```

**Cross-references:**

- Use `[AnnData-usage]` for user-facing documentation
- Link related functions: `@seealso [read_h5ad()], [write_h5ad()]`

### Vignettes

- `vignettes/anndataR.Rmd`: Main usage guide
- `vignettes/software_design.Rmd`: Architecture documentation
- `vignettes/usage_*.Rmd`: Format-specific tutorials
- `vignettes/known_issues.Rmd`: Documented limitations

## CI/CD Considerations

**Bioconductor requirements:**

- Must pass `R CMD check --as-cran` with no errors/warnings
- BiocCheck compliance
- All examples must run successfully
- Conditional package usage (see `check_requires` pattern)

## Quick Reference

### Creating AnnData Objects

```r
# In-memory (default)
adata <- AnnData(
  X = matrix(1:15, 3L, 5L),
  obs = data.frame(row.names = LETTERS[1:3]),
  var = data.frame(row.names = letters[1:5])
)

# HDF5-backed
adata <- HDF5AnnData$new("path.h5ad", mode = "w")
adata$X <- matrix(1:15, 3L, 5L)
# ... set other slots ...

# From file
adata <- read_h5ad("file.h5ad")
```

### Conversion Examples

```r
# SingleCellExperiment → AnnData
adata <- from_SingleCellExperiment(sce, X_name = "counts")

# AnnData → SingleCellExperiment
sce <- as_SingleCellExperiment(adata)

# Seurat → AnnData
adata <- from_Seurat(seurat_obj)

# AnnData → Seurat
seurat <- as_Seurat(adata)

# S4 coercion (if registered)
sce <- as(adata, "SingleCellExperiment")
```

### Accessing Data

```r
# Dimensions (using S3 methods)
dim(adata)         # [n_obs, n_vars]
nrow(adata)        # n_obs
ncol(adata)        # n_vars
dimnames(adata)    # list(obs_names, var_names)
rownames(adata)    # obs_names
colnames(adata)    # var_names

# Or using R6 methods
adata$n_obs()      # Number of observations
adata$n_vars()     # Number of variables

# Main matrix
adata$X            # Read
adata$X <- mat     # Write

# Metadata
adata$obs          # Observation metadata
adata$var          # Variable metadata
adata$uns          # Unstructured metadata

# Additional matrices
adata$layers[["raw"]]      # Named layer
adata$obsm[["X_pca"]]      # Observation matrix (PCA)
adata$obsp[["distances"]]  # Pairwise matrix
```

### Subsetting with AnnDataView

```r
# Subsetting returns AnnDataView (lazy, no data copied)
view <- adata[1:10, ]                     # Subset observations
view <- adata[, c("gene1", "gene2")]      # Subset variables
view <- adata[adata$obs$cell_type == "T", ]  # Conditional subsetting

# Convert to concrete implementation to apply changes
result <- view$as_InMemoryAnnData()
result <- view$as_HDF5AnnData("output.h5ad")
```

## Questions to Ask When Contributing

1. **Does this need Python compatibility?** → Add roundtrip test
2. **Requires optional package?** → Use `check_requires()` + conditional tests
3. **Modifying validation?** → Check AbstractAnnData validators
4. **New matrix type?** → Add to `generate_dataset()` for testing
5. **Changing uns handling?** → Consider Python dict requirements
6. **HDF5 changes?** → Verify with `h5diff` against Python output
7. **New conversion feature?** → Update both `from_*()` and `as_*()` paths

## Additional Resources

- [AnnData spec](https://anndata.readthedocs.io/): Official Python documentation
- `vignettes/software_design.Rmd`: Detailed architecture diagrams
- `inst/known_issues.yaml`: Tracked bugs and workarounds
- [Bioconductor submission guidelines](https://bioconductor.org/developers/)
