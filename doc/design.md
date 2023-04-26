# Design document

## Proposed interface

``` r
library(magicann)/library(annmu)

# read from h5ad/h5mu file
adata <- read_h5ad("dataset.h5ad")
adata <- read_h5ad("dataset.h5ad", backed = TRUE)
mdata <- read_h5mu("dataset.h5mu")
mdata <- read_h5mu("dataset.h5mu", backed = TRUE)

# optional feature 1: python-like interface
adata$X
adata$obs
adata$var

# optional feature 2: bioconductor-like interface
rowData(adata)
colData(adata)
reducedDimNames(adata)

# converters from/to sce
sce <- adata$to_sce()
from_sce(sce)

# optional feature 3: converters from/to seurat
seu <- adata$to_seurat()
from_seurat(seu)

# optional feature 4: converters from/to SOMA
som <- adata$to_soma()
from_soma(som)
```

## Class diagram

``` mermaid
classDiagram
  class AbstractAnnData {
    X: Matrix
    layers: List[Matrix]
    obs: DataFrame
    var: DataFrame
    obsp: List[Matrix]
    varp: List[Matrix]
    obsm: List[Matrix]
    varm: List[Matrix]
    uns: List
    n_obs: int
    n_vars: int
    obs_names: Array[String]
    var_names: Array[String]
    subset(...): AbstractAnnData
    to_sce(): SingleCellExperiment
    to_seurat(): Seurat
    write_h5ad(): Unit
  }
  AbstractAnnData <|-- BackedH5AD
  AbstractAnnData <|-- BackedZarr
  AbstractAnnData <|-- InMemoryAnnData

  class Package {
    read_h5ad(): AbstractAnnData
    read_h5mu(): AbstractMuData
  } 
```

## OO-framework

S4, RC, or R6?

- S4 offers formal class definitions and multiple dispatch, making it suitable for complex projects, but may be verbose and slower compared to other systems.
- RC provides reference semantics, familiar syntax, and encapsulation, yet it is less popular and can have performance issues. 
- R6 presents a simple and efficient OOP system with reference semantics and growing popularity, but lacks multiple dispatch and the formality of S4. 

Choosing an OOP system depends on the project requirements, developer familiarity, and desired balance between formality, performance, and ease of use.

## Approach

* Implement inheritance objects for `AbstractAnnData`, `BackedH5AD`, `InMemoryAnnData`
* Only containing `X`, `obs`, `var` for now
* `read_h5ad`, `write_h5ad`
* `$to_sce()` and `$to_seurat()`
* Explore adding `BackedZarr`