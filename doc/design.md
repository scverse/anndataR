# Design document

## Proposed interface

``` r
library(anndataR)

# read from h5ad/h5mu file
adata <- read_h5ad("dataset.h5ad", backend = "HDF5AnnData")
adata <- read_h5ad("dataset.h5ad", backend = "InMemoryAnnData")

# anndata-like interface (the Python package)
adata$X
adata$obs
adata$var

# optional feature 1: S3 helper functions for a base R-like interface
adata[1:10, 2:30]
dim(adata)
dimnames(adata)
as.matrix(adata, layer = NULL)
as.matrix(adata, layer = "counts")
t(adata)

# optional feature 2: S3 helper functions for a bioconductor-like interface
rowData(adata)
colData(adata)
reducedDimNames(adata)

# converters from/to sce
sce <- adata$to_SingleCellExperiment()
from_sce(sce)

# optional feature 3: converters from/to Seurat
seu <- adata$to_Seurat()
from_seurat(seu)

# optional feature 4: converters from/to SOMA
som <- adata$to_SOMA()
from_soma(som)
```

## Class diagram

``` mermaid
classDiagram
  class AbstractAnnData {
    *X: Matrix
    *layers: List[Matrix]
    *obs: DataFrame
    *var: DataFrame
    *obsp: List[Matrix]
    *varp: List[Matrix]
    *obsm: List[Matrix]
    *varm: List[Matrix]
    *uns: List
    *n_obs: int
    *n_vars: int
    *obs_names: Array[String]
    *var_names: Array[String]
    *subset(...): AbstractAnnData
    *write_h5ad(): Unit

    to_SingleCellExperiment(): SingleCellExperiment
    to_Seurat(): Seurat

    to_HDF5AnnData(): HDF5AnnData
    to_ZarrAnnData(): ZarrAnnData
    to_InMemoryAnnData(): InMemoryAnnData
  }

  AbstractAnnData <|-- HDF5AnnData
  class HDF5AnnData {
    init(h5file): HDF5AnnData
  }

  AbstractAnnData <|-- ZarrAnnData
  class ZarrAnnData {
    init(zarrFile): ZarrAnnData
  }

  AbstractAnnData <|-- InMemoryAnnData
  class InMemoryAnnData {
    init(X, obs, var, shape, ...): InMemoryAnnData
  }

  AbstractAnnData <|-- ReticulateAnnData
  class ReticulateAnnData {
    init(pyobj): ReticulateAnnData
  }

  class anndataR {
    read_h5ad(path, backend): Either[AbstractAnnData, SingleCellExperiment, Seurat]
  }
  anndataR --> AbstractAnnData
```

Notation:

- `X: Matrix` - variable `X` is of type `Matrix`
- `*X: Matrix` - variable `X` is abstract
- `to_SingleCellExperiment(): SingleCellExperiment` - function
  `to_SingleCellExperiment` returns object of type
  `SingleCellExperiment`
- `*to_SingleCellExperiment()` - function `to_SingleCellExperiment` is
  abstract

## OO-framework

S4, RC, or R6?

- S4 offers formal class definitions and multiple dispatch, making it
  suitable for complex projects, but may be verbose and slower compared
  to other systems.
- RC provides reference semantics, familiar syntax, and encapsulation,
  yet it is less popular and can have performance issues.
- R6 presents a simple and efficient OOP system with reference semantics
  and growing popularity, but lacks multiple dispatch and the formality
  of S4.

Choosing an OOP system depends on the project requirements, developer
familiarity, and desired balance between formality, performance, and
ease of use.

## Approach

- Implement inheritance objects for `AbstractAnnData`, `HDF5AnnData`,
  `InMemoryAnnData`
- Only containing `X`, `layers`, `obs`, `var` for now
- Implement base R S3 generics
- Implement `read_h5ad()`, `$write_h5ad()`
- Implement `$to_SingleCellExperiment()`
- Add simple unit tests

Optional:

- Add more fields (obsp, obsm, varp, varm, …) –\> see class diagram
- Start implementing MuData
- Implement `$to_Seurat()`
- Implement `ZarrAnnData`
- Implement `ReticulateAnnData`
- Implement Bioconductor S3 generics
