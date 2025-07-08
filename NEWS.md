# anndataR 0.1.0.9003

- Simplify & update vignetttes

# anndataR 0.99.0

## New features

* PR #158: Change package version from 0.0.0.9000 --> 0.99.0 to align with Bioc devel 
    versioning standards.
  - Update DESCRIPTION file to reflect current release R version used by Bio: 
    R (>= 3.2.0) --> R (>= 4.3.0).
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
    
  