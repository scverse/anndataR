# Proposed interface

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

## Alternative notation

- `ReadH5AD`, `ReadH5MU`, `WriteH5AD`, `WriteH5MU`, `adata$ToSCE()`,
  `adata$ToSeurat()`

## Class diagram

``` mermaid
classDiagram
  class AbstractAnnData {
    +Matrix X
    +List[Matrix] layers
    +DataFrame obs
    +DataFrame var
    +List[Matrix] obsp
    +List[Matrix] varp
    +List[Matrix] obsm
    +List[Matrix] varm
    +List uns
    +int n_obs
    +int n_vars
    +Array[String] obs_names
    +Array[String] var_names
    +subset(...): AbstractAnnData
    +toSCE(): SingleCellExperiment
    +toSeurat(): Seurat
  }
  AbstractAnnData &lt;|-- BackedH5AD
  AbstractAnnData &lt;|-- BackedZarr
  AbstractAnnData &lt;|-- InMemoryAnnData
```

## OO-framework

S4 or R6?
