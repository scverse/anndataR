# Generate a mock dataset

Generate a mock synthetic dataset with different types of columns and
layers. This is primarily designed for use in tests, examples, vignettes
and other documentation but is also provided to users for creating
reproducible examples.

Use `get_generator_types()` to get the available types for each slot.

## Usage

``` r
generate_dataset(
  n_obs = 10L,
  n_vars = 20L,
  x_type = "numeric_matrix",
  layer_types = get_generator_types(slot = "layers"),
  obs_types = get_generator_types(slot = "obs"),
  var_types = get_generator_types(slot = "var"),
  obsm_types = get_generator_types(slot = "obsm"),
  varm_types = get_generator_types(slot = "varm"),
  obsp_types = get_generator_types(slot = "obsp"),
  varp_types = get_generator_types(slot = "varp"),
  uns_types = get_generator_types(slot = "uns"),
  example = FALSE,
  format = c("list", "AnnData", "SingleCellExperiment", "Seurat")
)

get_generator_types(example = FALSE, slot = NULL)
```

## Arguments

- n_obs:

  Number of observations to generate

- n_vars:

  Number of variables to generate

- x_type:

  Type of matrix to generate for `X`

- layer_types:

  Types of matrices to generate for `layers`

- obs_types:

  Types of vectors to generate for `obs`

- var_types:

  Types of vectors to generate for `var`

- obsm_types:

  Types of matrices to generate for `obsm`

- varm_types:

  Types of matrices to generate for `varm`

- obsp_types:

  Types of matrices to generate for `obsp`

- varp_types:

  Types of matrices to generate for `varp`

- uns_types:

  Types of objects to generate for `uns`

- example:

  If `TRUE`, the types will be overridden with a small subset of types.
  This is useful for documentation.

- format:

  Object type to output, one of "list", "AnnData",
  "SingleCellExperiment", or "Seurat".

- slot:

  Which slot to return types for, if `NULL` a named list of all slots is
  returned

## Value

For `generate_dataset()`, an object as defined by `output` containing
the generated dataset

For `get_generator_types()`, a named list of character vectors or a
single character vector if `slot` is not `NULL`

## Details

To generate no data for a given slot, set the matching type argument to
`NULL` or an empty vector, e.g. `obs_types = c()` will generate an empty
`obs` data frame.

When generating `SingleCellExperiment` or `Seurat` objects, only some of
the generated slots will be included in the output object. To generate a
more complete object, use `format = "AnnData"` followed by
`adata$as_SingleCellExperiment()` or `adata$as_Seurat()`.

Use `get_generator_types()` to get a list of the available types for
each slot, or for a specific slot by setting `slot`. If
`example = TRUE`, only the example types are returned.

## Examples

``` r
# Generate all types as a list
dummy <- generate_dataset()

# Generate the example types
dummy_example <- generate_dataset(example = TRUE)

# Generate an AnnData
dummy_anndata <- generate_dataset(format = "AnnData", example = TRUE)

# Generate a SingleCellExperiment
if (rlang::is_installed("SingleCellExperiment")) {
  dummy_sce <- generate_dataset(format = "SingleCellExperiment", example = TRUE)
}

# Generate a Seurat object
if (rlang::is_installed("SeuratObject")) {
  dummy_seurat <- generate_dataset(format = "Seurat", example = TRUE)
}

# Get all available generator types
get_generator_types()
#> $X
#>  [1] "numeric_matrix"           "numeric_dense"           
#>  [3] "numeric_csparse"          "numeric_rsparse"         
#>  [5] "numeric_matrix_with_nas"  "numeric_dense_with_nas"  
#>  [7] "numeric_csparse_with_nas" "numeric_rsparse_with_nas"
#>  [9] "integer_matrix"           "integer_csparse"         
#> [11] "integer_rsparse"          "integer_matrix_with_nas" 
#> [13] "integer_csparse_with_nas" "integer_rsparse_with_nas"
#> 
#> $layers
#>  [1] "numeric_matrix"           "numeric_dense"           
#>  [3] "numeric_csparse"          "numeric_rsparse"         
#>  [5] "numeric_matrix_with_nas"  "numeric_dense_with_nas"  
#>  [7] "numeric_csparse_with_nas" "numeric_rsparse_with_nas"
#>  [9] "integer_matrix"           "integer_csparse"         
#> [11] "integer_rsparse"          "integer_matrix_with_nas" 
#> [13] "integer_csparse_with_nas" "integer_rsparse_with_nas"
#> 
#> $obs
#>  [1] "character"               "integer"                
#>  [3] "factor"                  "factor_ordered"         
#>  [5] "logical"                 "numeric"                
#>  [7] "character_with_nas"      "integer_with_nas"       
#>  [9] "factor_with_nas"         "factor_ordered_with_nas"
#> [11] "logical_with_nas"        "numeric_with_nas"       
#> 
#> $var
#>  [1] "character"               "integer"                
#>  [3] "factor"                  "factor_ordered"         
#>  [5] "logical"                 "numeric"                
#>  [7] "character_with_nas"      "integer_with_nas"       
#>  [9] "factor_with_nas"         "factor_ordered_with_nas"
#> [11] "logical_with_nas"        "numeric_with_nas"       
#> 
#> $obsm
#>  [1] "numeric_matrix"           "numeric_dense"           
#>  [3] "numeric_csparse"          "numeric_rsparse"         
#>  [5] "numeric_matrix_with_nas"  "numeric_dense_with_nas"  
#>  [7] "numeric_csparse_with_nas" "numeric_rsparse_with_nas"
#>  [9] "integer_matrix"           "integer_csparse"         
#> [11] "integer_rsparse"          "integer_matrix_with_nas" 
#> [13] "integer_csparse_with_nas" "integer_rsparse_with_nas"
#> [15] "character"                "integer"                 
#> [17] "factor"                   "factor_ordered"          
#> [19] "logical"                  "numeric"                 
#> [21] "character_with_nas"       "integer_with_nas"        
#> [23] "factor_with_nas"          "factor_ordered_with_nas" 
#> [25] "logical_with_nas"         "numeric_with_nas"        
#> 
#> $varm
#>  [1] "numeric_matrix"           "numeric_dense"           
#>  [3] "numeric_csparse"          "numeric_rsparse"         
#>  [5] "numeric_matrix_with_nas"  "numeric_dense_with_nas"  
#>  [7] "numeric_csparse_with_nas" "numeric_rsparse_with_nas"
#>  [9] "integer_matrix"           "integer_csparse"         
#> [11] "integer_rsparse"          "integer_matrix_with_nas" 
#> [13] "integer_csparse_with_nas" "integer_rsparse_with_nas"
#> [15] "character"                "integer"                 
#> [17] "factor"                   "factor_ordered"          
#> [19] "logical"                  "numeric"                 
#> [21] "character_with_nas"       "integer_with_nas"        
#> [23] "factor_with_nas"          "factor_ordered_with_nas" 
#> [25] "logical_with_nas"         "numeric_with_nas"        
#> 
#> $obsp
#>  [1] "numeric_matrix"           "numeric_dense"           
#>  [3] "numeric_csparse"          "numeric_rsparse"         
#>  [5] "numeric_matrix_with_nas"  "numeric_dense_with_nas"  
#>  [7] "numeric_csparse_with_nas" "numeric_rsparse_with_nas"
#>  [9] "integer_matrix"           "integer_csparse"         
#> [11] "integer_rsparse"          "integer_matrix_with_nas" 
#> [13] "integer_csparse_with_nas" "integer_rsparse_with_nas"
#> 
#> $varp
#>  [1] "numeric_matrix"           "numeric_dense"           
#>  [3] "numeric_csparse"          "numeric_rsparse"         
#>  [5] "numeric_matrix_with_nas"  "numeric_dense_with_nas"  
#>  [7] "numeric_csparse_with_nas" "numeric_rsparse_with_nas"
#>  [9] "integer_matrix"           "integer_csparse"         
#> [11] "integer_rsparse"          "integer_matrix_with_nas" 
#> [13] "integer_csparse_with_nas" "integer_rsparse_with_nas"
#> 
#> $uns
#>  [1] "scalar_character"               "scalar_integer"                
#>  [3] "scalar_factor"                  "scalar_factor_ordered"         
#>  [5] "scalar_logical"                 "scalar_numeric"                
#>  [7] "scalar_character_with_nas"      "scalar_integer_with_nas"       
#>  [9] "scalar_factor_with_nas"         "scalar_factor_ordered_with_nas"
#> [11] "scalar_logical_with_nas"        "scalar_numeric_with_nas"       
#> [13] "vec_character"                  "vec_integer"                   
#> [15] "vec_factor"                     "vec_factor_ordered"            
#> [17] "vec_logical"                    "vec_numeric"                   
#> [19] "vec_character_with_nas"         "vec_integer_with_nas"          
#> [21] "vec_factor_with_nas"            "vec_factor_ordered_with_nas"   
#> [23] "vec_logical_with_nas"           "vec_numeric_with_nas"          
#> [25] "df_character"                   "df_integer"                    
#> [27] "df_factor"                      "df_factor_ordered"             
#> [29] "df_logical"                     "df_numeric"                    
#> [31] "df_character_with_nas"          "df_integer_with_nas"           
#> [33] "df_factor_with_nas"             "df_factor_ordered_with_nas"    
#> [35] "df_logical_with_nas"            "df_numeric_with_nas"           
#> [37] "mat_numeric_matrix"             "mat_numeric_dense"             
#> [39] "mat_numeric_csparse"            "mat_numeric_rsparse"           
#> [41] "mat_numeric_matrix_with_nas"    "mat_numeric_dense_with_nas"    
#> [43] "mat_numeric_csparse_with_nas"   "mat_numeric_rsparse_with_nas"  
#> [45] "mat_integer_matrix"             "mat_integer_csparse"           
#> [47] "mat_integer_rsparse"            "mat_integer_matrix_with_nas"   
#> [49] "mat_integer_csparse_with_nas"   "mat_integer_rsparse_with_nas"  
#> [51] "list"                          
#> 

# Get generator types for a specific slot
get_generator_types(slot = "obs")
#>  [1] "character"               "integer"                
#>  [3] "factor"                  "factor_ordered"         
#>  [5] "logical"                 "numeric"                
#>  [7] "character_with_nas"      "integer_with_nas"       
#>  [9] "factor_with_nas"         "factor_ordered_with_nas"
#> [11] "logical_with_nas"        "numeric_with_nas"       

# Get generator types used when example = TRUE
get_generator_types(example = TRUE)
#> $X
#> [1] "numeric_matrix"
#> 
#> $layers
#> [1] "numeric_matrix"  "numeric_dense"   "numeric_csparse"
#> 
#> $obs
#> [1] "character" "integer"   "factor"   
#> 
#> $var
#> [1] "character" "integer"   "factor"   
#> 
#> $obsm
#> [1] "numeric_matrix"  "numeric_dense"   "numeric_csparse"
#> 
#> $varm
#> [1] "numeric_matrix"  "numeric_dense"   "numeric_csparse"
#> 
#> $obsp
#> [1] "numeric_matrix"  "numeric_dense"   "numeric_csparse"
#> 
#> $varp
#> [1] "numeric_matrix"  "numeric_dense"   "numeric_csparse"
#> 
#> $uns
#> [1] "scalar_character"   "vec_integer"        "df_logical"        
#> [4] "mat_numeric_matrix"
#> 
```
