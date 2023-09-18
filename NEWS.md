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