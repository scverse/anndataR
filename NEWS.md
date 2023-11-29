<<<<<<< HEAD
# anndataR 0.99.1

## New features

* Harmonise all references to `h5ad_file` (`file`,`h5ad_path`,`ad`)
* Harmonise all references to `ad` (`adata`)
* Use `@returns` in Roxygen notes to allow multi-line notes.
* Export helper function `dummy_data`
  - Makes examples much less verbose and more consistent.
  - Include "obsm","varm","obsp","varp" and "layers" wherever possible.
* *Seurat.R*
  - Add/extract `DimRed` elements whenever available.
  - Add unit tests for `from_Seurat`.
* New functions
  - `graphs_to_anndata`:  When converting Seurat --> AnnData, 
    store `obsp` and `varp` anndata slots from `graphs`. Does some inference 
    about whether each matrix should go into the `obsp` vs. `varp` slot based 
    on whether they have the same dimensions as the number of rows. 
  - `anndata_to_graphs`: When converting AnnData --> Seurat, 
    store `obsp` and `varp` anndata slots in `graphs`. 
* *DimReduc.R*
  - New functions for converting Seurat's DimReduc objects.
  

## Bug fixes

* *DESCRIPTION*:
  - Add "email" tag to `Authors`.
  - Add `BiocViews`
* Fix `BiocCheck` WARNINGS:
  - `Empty or missing \value sections found in man pages.`
* Fix `BiocCheck` NOTES:
  - Make all vignettes use `BiocStyle`
  - `Update R version dependency from 3.2.0 to 4.3.0.`
  - `'sessionInfo' not found in vignette(s)`
  - `Avoid the use of 'paste' in condition signals` 
  - `Avoid 'cat' and 'print' outside of 'show' methods`
  - `Consider adding runnable examples to man pages that document exported objects.`
  - Added helper function `messager` to avoid notes about pasting inside messages.
    Can also easily turn off verbosity by setting the env var "ANNDATAR_VERBOSE".
* Fix all unit tests given updates.

# anndataR 0.99.0

## New features

- Add `rworkflows` CI.
- Update *lint.yaml* to use `actions/checkout@v4` (which has less issues).
- New function `setup_conda` automatically installs miniconda 
  and sets up conda env: #97
- Change version to Bioc-recommended devel version: 0.99.0

# anndataR 0.1.0

## New features

Initial release of anndataR, providing support for working with AnnData objects in R.

Feature list:

* Slots:
  - X
  - layers
  - obs
  - obs_names
  - var
  - var_names

* Backends:
  - HDF5AnnData
  - InMemoryAnnData

* Converters:
  - SingleCellExperiment
  - Seurat
