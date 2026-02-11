# S3 Methods for AbstractAnnData Objects

These S3 methods provide standard R interfaces for AbstractAnnData
objects, making them behave like native R objects with familiar syntax.

## Usage

``` r
# S3 method for class 'AbstractAnnData'
dim(x)

# S3 method for class 'AbstractAnnData'
nrow(x)

# S3 method for class 'AbstractAnnData'
ncol(x)

# S3 method for class 'AbstractAnnData'
dimnames(x)

# S3 method for class 'AbstractAnnData'
dimnames(x) <- value

# S3 method for class 'AbstractAnnData'
x[i, j, drop = TRUE, ...]
```

## Arguments

- x:

  An AbstractAnnData object

- value:

  For `dimnames<-`: A list of two character vectors (obs_names,
  var_names)

- i:

  Row indices (observations). Can be numeric, logical, or character.

- j:

  Column indices (variables). Can be numeric, logical, or character.

- drop:

  Ignored (for compatibility with generic)

- ...:

  Additional arguments passed to methods

## Value

- `dim`: Numeric vector of length 2 (n_obs, n_vars)

- `nrow`, `ncol`: Integer count

- `dimnames`: List with obs_names and var_names

- `rownames`, `colnames`: Character vector

- `dimnames<-`, `rownames<-`, `colnames<-`: The modified object
  (invisibly)

- `[`: A AnnDataView object with the specified subset

## Details

**Subsetting behaviour**: The `[` method supports logical, integer, and
character subsetting for both observations (rows) and variables
(columns). However, unlike standard R behaviour:

- Logical vectors are **not recycled** and must have the exact same
  length as the dimension being subset

- **Negative indices are not supported** (R's "exclude these" syntax)

These design choices ensure clear and predictable subsetting behaviour
for biological data matrices, avoiding potential confusion from
accidental recycling or exclusion patterns.

The following S3 methods are available:

- `dim(x)`: Get dimensions (n_obs, n_vars), equivalent to `x$shape()`

- `nrow(x)`: Get number of observations, equivalent to `x$n_obs()`

- `ncol(x)`: Get number of variables, equivalent to `x$n_vars()`

- `dimnames(x)`: Get dimension names, (obs_names, var_names)

- `rownames(x)`: Get observation names, equivalent to `x$obs_names()`

- `colnames(x)`: Get variable names, equivalent to `x$var_names()`

- `dimnames(x) <- value`: Set dimension names

- `rownames(x) <- value`: Set observation names, equivalent to
  `x$obs_names() <- ...`

- `colnames(x) <- value`: Set variable names, equivalent to
  `x$var_names() <- ...`

- `x[i, j]`: Subset observations and/or variables

## Examples

``` r
# Create example data
ad <- generate_dataset(n_obs = 100, n_vars = 50, format = "AnnData")

# Standard R methods work
dim(ad)
#> [1] 100  50
nrow(ad)
#> [1] 100
ncol(ad)
#> [1] 50
dimnames(ad)
#> [[1]]
#>   [1] "cell1"   "cell2"   "cell3"   "cell4"   "cell5"   "cell6"   "cell7"  
#>   [8] "cell8"   "cell9"   "cell10"  "cell11"  "cell12"  "cell13"  "cell14" 
#>  [15] "cell15"  "cell16"  "cell17"  "cell18"  "cell19"  "cell20"  "cell21" 
#>  [22] "cell22"  "cell23"  "cell24"  "cell25"  "cell26"  "cell27"  "cell28" 
#>  [29] "cell29"  "cell30"  "cell31"  "cell32"  "cell33"  "cell34"  "cell35" 
#>  [36] "cell36"  "cell37"  "cell38"  "cell39"  "cell40"  "cell41"  "cell42" 
#>  [43] "cell43"  "cell44"  "cell45"  "cell46"  "cell47"  "cell48"  "cell49" 
#>  [50] "cell50"  "cell51"  "cell52"  "cell53"  "cell54"  "cell55"  "cell56" 
#>  [57] "cell57"  "cell58"  "cell59"  "cell60"  "cell61"  "cell62"  "cell63" 
#>  [64] "cell64"  "cell65"  "cell66"  "cell67"  "cell68"  "cell69"  "cell70" 
#>  [71] "cell71"  "cell72"  "cell73"  "cell74"  "cell75"  "cell76"  "cell77" 
#>  [78] "cell78"  "cell79"  "cell80"  "cell81"  "cell82"  "cell83"  "cell84" 
#>  [85] "cell85"  "cell86"  "cell87"  "cell88"  "cell89"  "cell90"  "cell91" 
#>  [92] "cell92"  "cell93"  "cell94"  "cell95"  "cell96"  "cell97"  "cell98" 
#>  [99] "cell99"  "cell100"
#> 
#> [[2]]
#>  [1] "gene1"  "gene2"  "gene3"  "gene4"  "gene5"  "gene6"  "gene7"  "gene8" 
#>  [9] "gene9"  "gene10" "gene11" "gene12" "gene13" "gene14" "gene15" "gene16"
#> [17] "gene17" "gene18" "gene19" "gene20" "gene21" "gene22" "gene23" "gene24"
#> [25] "gene25" "gene26" "gene27" "gene28" "gene29" "gene30" "gene31" "gene32"
#> [33] "gene33" "gene34" "gene35" "gene36" "gene37" "gene38" "gene39" "gene40"
#> [41] "gene41" "gene42" "gene43" "gene44" "gene45" "gene46" "gene47" "gene48"
#> [49] "gene49" "gene50"
#> 
rownames(ad)
#>   [1] "cell1"   "cell2"   "cell3"   "cell4"   "cell5"   "cell6"   "cell7"  
#>   [8] "cell8"   "cell9"   "cell10"  "cell11"  "cell12"  "cell13"  "cell14" 
#>  [15] "cell15"  "cell16"  "cell17"  "cell18"  "cell19"  "cell20"  "cell21" 
#>  [22] "cell22"  "cell23"  "cell24"  "cell25"  "cell26"  "cell27"  "cell28" 
#>  [29] "cell29"  "cell30"  "cell31"  "cell32"  "cell33"  "cell34"  "cell35" 
#>  [36] "cell36"  "cell37"  "cell38"  "cell39"  "cell40"  "cell41"  "cell42" 
#>  [43] "cell43"  "cell44"  "cell45"  "cell46"  "cell47"  "cell48"  "cell49" 
#>  [50] "cell50"  "cell51"  "cell52"  "cell53"  "cell54"  "cell55"  "cell56" 
#>  [57] "cell57"  "cell58"  "cell59"  "cell60"  "cell61"  "cell62"  "cell63" 
#>  [64] "cell64"  "cell65"  "cell66"  "cell67"  "cell68"  "cell69"  "cell70" 
#>  [71] "cell71"  "cell72"  "cell73"  "cell74"  "cell75"  "cell76"  "cell77" 
#>  [78] "cell78"  "cell79"  "cell80"  "cell81"  "cell82"  "cell83"  "cell84" 
#>  [85] "cell85"  "cell86"  "cell87"  "cell88"  "cell89"  "cell90"  "cell91" 
#>  [92] "cell92"  "cell93"  "cell94"  "cell95"  "cell96"  "cell97"  "cell98" 
#>  [99] "cell99"  "cell100"
colnames(ad)
#>  [1] "gene1"  "gene2"  "gene3"  "gene4"  "gene5"  "gene6"  "gene7"  "gene8" 
#>  [9] "gene9"  "gene10" "gene11" "gene12" "gene13" "gene14" "gene15" "gene16"
#> [17] "gene17" "gene18" "gene19" "gene20" "gene21" "gene22" "gene23" "gene24"
#> [25] "gene25" "gene26" "gene27" "gene28" "gene29" "gene30" "gene31" "gene32"
#> [33] "gene33" "gene34" "gene35" "gene36" "gene37" "gene38" "gene39" "gene40"
#> [41] "gene41" "gene42" "gene43" "gene44" "gene45" "gene46" "gene47" "gene48"
#> [49] "gene49" "gene50"

# Set names using dimnames
dimnames(ad) <- list(
  paste0("cell_", 1:nrow(ad)),
  paste0("gene_", 1:ncol(ad))
)

# Or set names individually (uses dimnames<- internally)
rownames(ad) <- paste0("cell_", 1:nrow(ad))
colnames(ad) <- paste0("gene_", 1:ncol(ad))

# Subsetting creates AnnDataView
subset_ad <- ad[1:10, 1:5]
subset_ad <- ad[rep(c(TRUE, FALSE), length.out = nrow(ad)), ]  # logical subsetting (no recycling)
subset_ad <- ad[c("cell_1", "cell_2"), c("gene_1", "gene_2")]  # name subsetting
```
