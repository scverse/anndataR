#' Read an element
#'
#' Read a single element from an H5AD file or Zarr store
#'
#' @param path Path to the H5AD file or Zarr store to read from
#' @param name Path of the element inside the file or store, e.g. `"obs"` or
#'   `"obsm/X_pca"`
#' @param backed Whether to return array and matrix elements as
#'   [DelayedArray::DelayedArray] objects instead of reading them into memory.
#'   Only supported for H5AD files.
#' @param format The format of `path`. One of:
#'
#'   * `"auto"`: Detect the format from `path`. A file that is an HDF5 file is
#'     read as H5AD and a directory with Zarr metadata is read as Zarr.
#'   * `"h5ad"`: Read `path` as an H5AD file
#'   * `"zarr"`: Read `path` as a Zarr store
#' @param ... Not currently used. Unused arguments are an error.
#'
#' @return The element, as the matching R type
#'
#' @export
#'
#' @family element readers and writers
#'
#' @examples
#' h5ad_file <- system.file("extdata", "example.h5ad", package = "anndataR")
#'
#' if (requireNamespace("rhdf5", quietly = TRUE)) {
#'   # Read the obs data frame
#'   obs <- read_elem(h5ad_file, "obs")
#'
#'   # Read a single obsm matrix
#'   pca <- read_elem(h5ad_file, "obsm/X_pca")
#' }
read_elem <- function(
  path,
  name,
  backed = FALSE,
  format = c("auto", "h5ad", "zarr"),
  ...
) {
  rlang::check_dots_used()

  check_elem_name(name)
  if (!rlang::is_bool(backed)) {
    cli_abort("{.arg backed} must be {.val {TRUE}} or {.val {FALSE}}")
  }
  if (rlang::is_string(path) && !file.exists(path)) {
    cli_abort("{.path {path}} does not exist")
  }

  format <- detect_elem_format(path, format)

  switch(
    format,
    h5ad = {
      check_requires("Reading H5AD elements", "rhdf5", where = "Bioc")
      check_requires("Reading H5AD elements", "withr", where = "CRAN")
      if (backed) {
        check_requires(
          "Reading backed H5AD elements",
          c("DelayedArray", "HDF5Array"),
          where = "Bioc"
        )
      }
      .read_elem_h5ad(HDF5File$new(path), name, backed = backed, ...)
    },
    zarr = {
      check_requires("Reading Zarr elements", "Rarr", where = "Bioc")
      if (backed) {
        # TODO: Support backed Zarr elements once #507 is merged
        cli_abort("{.code backed = TRUE} is not yet supported for Zarr stores")
      }
      .read_elem_zarr(path, name, ...)
    }
  )
}

#' Read an element from an H5AD file
#'
#' @param hdf5_file An `HDF5File` object
#' @param name Path of the element inside the file
#' @param backed Whether to return a `DelayedArray` for array/matrix elements
#' @param ... Passed to `read_h5ad_element()`
#' @param call The calling environment, used in error messages
#'
#' @return The element
#'
#' @noRd
.read_elem_h5ad <- function(
  hdf5_file,
  name,
  backed = FALSE,
  ...,
  call = rlang::caller_env()
) {
  name <- hdf5_root_path("/", name)

  # H5Lexists() can't check the root group, which always exists
  if (name != "/" && !hdf5_path_exists(hdf5_file, name)) {
    cli_abort(
      "Element {.val {name}} not found in {.path {hdf5_file$path}}",
      call = call
    )
  }

  read_h5ad_element(
    hdf5_file,
    name,
    stop_on_error = TRUE,
    backed = backed,
    ...
  )
}

#' Read an element from a Zarr store
#'
#' @param store Path to the Zarr store
#' @param name Path of the element inside the store
#' @param ... Passed to `read_zarr_element()`
#' @param call The calling environment, used in error messages
#'
#' @return The element
#'
#' @noRd
.read_elem_zarr <- function(store, name, ..., call = rlang::caller_env()) {
  name <- gsub("^/+|/+$", "", name)
  if (name == "") {
    name <- "/"
  }

  if (!zarr_path_exists(store, name)) {
    cli_abort("Element {.val {name}} not found in {.path {store}}", call = call)
  }

  read_zarr_element(
    store,
    name,
    stop_on_error = TRUE,
    ...
  )
}

#' Check an element name
#'
#' Check that `name` is a valid element name for `read_elem()`/`write_elem()`
#'
#' @param name The element name to check
#' @param call The calling environment, used in error messages
#'
#' @return `name` (invisibly)
#'
#' @noRd
check_elem_name <- function(name, call = rlang::caller_env()) {
  if (!rlang::is_string(name) || is.na(name)) {
    cli_abort(
      "{.arg name} must be a single string, not {.obj_type_friendly {name}}",
      call = call
    )
  }

  invisible(name)
}
