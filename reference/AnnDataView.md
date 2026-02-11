# A View of an AnnData Object

A lazy view of an AnnData object that allows applying subsetting
operations without immediately executing them. Subsetting is applied
when converting to a concrete AnnData implementation (InMemoryAnnData,
HDF5AnnData) or other formats (SingleCellExperiment, Seurat).

## Value

A `AnnDataView` object

## See also

[AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)
for details on creating and using `AnnData` objects

Other AnnData classes:
[`AbstractAnnData`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.md),
[`HDF5AnnData`](https://anndataR.data-intuitive.com/reference/HDF5AnnData.md),
[`InMemoryAnnData`](https://anndataR.data-intuitive.com/reference/InMemoryAnnData.md),
[`ReticulateAnnData`](https://anndataR.data-intuitive.com/reference/ReticulateAnnData.md)

## Super class

[`anndataR::AbstractAnnData`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.md)
-\> `AnnDataView`

## Active bindings

- `X`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `layers`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `obs`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `var`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `obs_names`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `var_names`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `obsm`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `varm`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `obsp`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `varp`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

- `uns`:

  See
  [AnnData-usage](https://anndataR.data-intuitive.com/reference/AnnData-usage.md)

## Methods

### Public methods

- [`AnnDataView$new()`](#method-AnnDataView-new)

- [`AnnDataView$subset()`](#method-AnnDataView-subset)

- [`AnnDataView$clone()`](#method-AnnDataView-clone)

Inherited methods

- [`anndataR::AbstractAnnData$as_HDF5AnnData()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-as_HDF5AnnData)
- [`anndataR::AbstractAnnData$as_InMemoryAnnData()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-as_InMemoryAnnData)
- [`anndataR::AbstractAnnData$as_ReticulateAnnData()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-as_ReticulateAnnData)
- [`anndataR::AbstractAnnData$as_Seurat()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-as_Seurat)
- [`anndataR::AbstractAnnData$as_SingleCellExperiment()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-as_SingleCellExperiment)
- [`anndataR::AbstractAnnData$layers_keys()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-layers_keys)
- [`anndataR::AbstractAnnData$n_obs()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-n_obs)
- [`anndataR::AbstractAnnData$n_vars()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-n_vars)
- [`anndataR::AbstractAnnData$obs_keys()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-obs_keys)
- [`anndataR::AbstractAnnData$obsm_keys()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-obsm_keys)
- [`anndataR::AbstractAnnData$obsp_keys()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-obsp_keys)
- [`anndataR::AbstractAnnData$print()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-print)
- [`anndataR::AbstractAnnData$shape()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-shape)
- [`anndataR::AbstractAnnData$uns_keys()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-uns_keys)
- [`anndataR::AbstractAnnData$var_keys()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-var_keys)
- [`anndataR::AbstractAnnData$varm_keys()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-varm_keys)
- [`anndataR::AbstractAnnData$varp_keys()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-varp_keys)
- [`anndataR::AbstractAnnData$write_h5ad()`](https://anndataR.data-intuitive.com/reference/AbstractAnnData.html#method-write_h5ad)

------------------------------------------------------------------------

### Method `new()`

Create a new AnnDataView object

#### Usage

    AnnDataView$new(base_adata, i, j)

#### Arguments

- `base_adata`:

  An existing AnnData object to create a view of

- `i`:

  Optional initial obs subset (logical, integer, or character vector)

- `j`:

  Optional initial var subset (logical, integer, or character vector)

#### Returns

A new `AnnDataView` object Subset the AnnDataView

------------------------------------------------------------------------

### Method [`subset()`](https://rdrr.io/r/base/subset.html)

#### Usage

    AnnDataView$subset(i, j)

#### Arguments

- `i`:

  Row indices (observations). Can be numeric, logical, or character.

- `j`:

  Column indices (variables). Can be numeric, logical, or character.

#### Returns

A new `AnnDataView` object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    AnnDataView$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# Create a base AnnData object
ad <- AnnData(
  X = matrix(1:15, 3L, 5L),
  obs = data.frame(row.names = LETTERS[1:3], cell_type = c("A", "B", "A")),
  var = data.frame(row.names = letters[1:5], gene_type = c("X", "Y", "X", "Y", "Z"))
)

# Create a view with lazy subsetting using S3 [ method
view <- ad[ad$obs$cell_type == "A", ad$var$gene_type %in% c("X", "Y")]

# Apply subsetting by converting to a concrete implementation
result <- view$as_InMemoryAnnData()
```
