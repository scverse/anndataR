#' @title Lazy Named List
#'
#' @description A lazy named list that loads elements on-demand to avoid
#' materializing large objects unnecessarily. Used internally for efficient
#' access to layers, obsm, varm, obsp, varp, uns slots.
#'
#' @keywords internal
AnnDataSlotList <- R6::R6Class(
  "AnnDataSlotList",
  public = list(
    get_keys_fn = NULL,
    set_keys_fn = NULL,
    get_value_fn = NULL,
    set_value_fn = NULL,
    set_values_fn = NULL,
    get_rownames_fn = NULL,
    get_colnames_fn = NULL,

    #' @description Create a new AnnDataSlotList
    #' @param get_keys_fn Function that returns all available keys: function() -> list of strings
    #' @param set_keys_fn Function to set all keys: function(keys) -> invisible()
    #' @param get_value_fn Function to get element by key: function(key) -> object
    #' @param set_value_fn Function to set element by key: function(key, value) -> invisible()
    #' @param set_values_fn Function to set multiple elements: function(named_list) -> invisible()
    #' @param get_rownames_fn An optional function to get the rownames the values should be aligned to: function() -> list of strings
    #' @param get_colnames_fn An optional function to get the colnames the values should be aligned to: function() -> list of strings
    initialize = function(
      get_keys_fn,
      set_keys_fn,
      get_value_fn,
      set_value_fn,
      set_values_fn,
      get_rownames_fn = NULL,
      get_colnames_fn = NULL
    ) {
      self$get_keys_fn <- get_keys_fn
      self$set_keys_fn <- set_keys_fn
      self$get_value_fn <- get_value_fn
      self$set_value_fn <- set_value_fn
      self$set_values_fn <- set_values_fn
      self$get_rownames_fn <- get_rownames_fn
      self$get_colnames_fn <- get_colnames_fn
    }
  )
)

