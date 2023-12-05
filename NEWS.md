# anndataR 0.99.0

## New features

* PR #159: Add dedicated internal functions for convert to/from Seurat DimReduc object.
  - Fix unit tests for Seurat converter.
  - Add missing extra final line to DESCRIPTION.
  - Make naming conventions consistent for internal functions: 
    - `.toseurat_check_obsvar_names` --> `.to_Seurat_check_obsvar_names`
  - Add helper internal function: `.to_Seurat_check_layer_names`

* PR #158: Change package version from 0.0.0.9000 --> 0.99.0 to align with Bioc devel 
    versioning standards.
  - Update DESCRIPTION file to reflect current release R version used by Bio: 
    R (>= 3.2.0) --> R (>= 4.0.0).
  - Reformat NEWS file to follow some conventions.

* Various PRs: Initial release of anndataR, providing support for working with
  AnnData objects in R. Feature list:
  - Slots:
    - X
    - layers
    - obs
    - obs_names
    - var
    - var_names
  - Backends:
    - HDF5AnnData
    - InMemoryAnnData
  - Converters:
    - SingleCellExperiment
    - Seurat
    
  