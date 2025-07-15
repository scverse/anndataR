convert_categorical <- function(categorical) {
  categorical <- reticulate::py_to_r(categorical)
  if (!inherits(categorical, "python.builtin.object")) {
    return(categorical)
  }

  categories <- reticulate::py_to_r(categorical$categories)
  codes <- reticulate::py_to_r(categorical$codes)
  ordered <- reticulate::py_to_r(categorical$ordered)
  is_na <- codes == -1L
  codes[is_na] <- 0L
  py_value <- factor(
    categories[codes + 1],
    levels = categories,
    ordered = ordered
  )
  py_value[is_na] <- NA

  py_value
}

convert_nullable_integer_array <- function(nullable_array) {
  nullable_array <- reticulate::py_to_r(nullable_array)
  if (!inherits(nullable_array, "python.builtin.object")) {
    return(nullable_array)
  }

  mask <- reticulate::py_to_r(reticulate::py_get_attr(nullable_array, "_mask"))
  data <- reticulate::py_to_r(reticulate::py_get_attr(nullable_array, "_data"))
  py_value <- as.integer(data)
  py_value[mask] <- NA

  py_value
}

convert_nullable_boolean_array <- function(nullable_array) {
  nullable_array <- reticulate::py_to_r(nullable_array)
  if (!inherits(nullable_array, "python.builtin.object")) {
    return(nullable_array)
  }

  mask <- reticulate::py_to_r(reticulate::py_get_attr(nullable_array, "_mask"))
  data <- reticulate::py_to_r(reticulate::py_get_attr(nullable_array, "_data"))
  py_value <- as.logical(data)
  py_value[mask] <- NA

  py_value
}
