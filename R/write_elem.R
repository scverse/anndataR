#' Write an element
#'
#' Write a single element to an H5AD file or Zarr store
#'
#' @param path Path to the H5AD file or Zarr store to write to. Created if it
#'   doesn't exist.
#' @param name Path of the element inside the file or store.
#' @param value The value to write
#' @param overwrite Whether to replace an element that already exists at
#'   `name`. If `FALSE` (default), an existing element is an error.
#' @param format The format of `path`. One of:
#'
#'   * `"auto"`: Detect the format from `path`. An existing HDF5 file is
#'     written as H5AD and an existing directory with Zarr metadata as Zarr. A
#'     new path is detected from its extension.
#'   * `"h5ad"`: Write `path` as an H5AD file
#'   * `"zarr"`: Write `path` as a Zarr store
#' @param ... Options for writing the element:
#'
#'   * `compression`: The compression to use. See [write_h5ad()] and
#'     [write_zarr()].
#'   * `chunk_size`: The target chunk size. See [write_h5ad()].
#'
#' @return `path` (invisibly)
#' @export
#'
#' @family element readers and writers
#'
#' @examples
#' h5ad_file <- tempfile(fileext = ".h5ad")
#'
#' if (requireNamespace("rhdf5", quietly = TRUE)) {
#'   # Write a data frame to a new file
#'   write_elem(h5ad_file, "uns/df", data.frame(a = 1:3, b = letters[1:3]))
#'
#'   # Replace it with a vector
#'   write_elem(h5ad_file, "uns/df", 1:10, overwrite = TRUE)
#'
#'   read_elem(h5ad_file, "uns/df")
#' }
write_elem <- function(
  path,
  name,
  value,
  overwrite = FALSE,
  format = c("auto", "h5ad", "zarr"),
  ...
) {
  check_elem_name(name)
  if (!rlang::is_bool(overwrite)) {
    cli_abort("{.arg overwrite} must be {.val {TRUE}} or {.val {FALSE}}")
  }

  format <- detect_elem_format(path, format)

  switch(
    format,
    h5ad = {
      check_requires("Writing H5AD elements", "rhdf5", where = "Bioc")
      check_requires("Writing H5AD elements", "withr", where = "CRAN")
      check_elem_dots(c("compression", "chunk_size"), ...)

      if (!file.exists(path)) {
        rhdf5::h5createFile(path)
      }

      .write_elem_h5ad(HDF5File$new(path), name, value, overwrite, ...)
    },
    zarr = {
      check_requires("Writing Zarr elements", "Rarr", where = "Bioc")
      check_elem_dots("compression", ...)

      path <- sub("/+$", "", path)
      zarr_format <- if (dir.exists(path)) {
        get_zarr_format(path)
      } else {
        zarr_format <- check_zarr_format(
          getOption("anndataR.zarr_format", 3L)
        )
        Rarr::write_zarr_group(path, "", zarr_version = zarr_format)
        zarr_format
      }

      .write_elem_zarr(path, name, value, overwrite, zarr_format, ...)
    }
  )

  invisible(path)
}

#' Write an element to an H5AD file
#'
#' @param hdf5_file An `HDF5File` object
#' @param name Path of the element inside the file
#' @param value The value to write
#' @param overwrite Whether to replace an existing element
#' @param ... Passed to `write_h5ad_element()`
#' @param call The calling environment, used in error messages
#'
#' @return `NULL` (invisibly)
#'
#' @noRd
.write_elem_h5ad <- function(
  hdf5_file,
  name,
  value,
  overwrite = FALSE,
  ...,
  call = rlang::caller_env()
) {
  name <- hdf5_root_path("/", name)

  # TODO: Allow writing AnnData objects to the root
  if (name == "/") {
    cli_abort("Writing to the root group is not supported yet", call = call)
  }

  if (hdf5_path_exists(hdf5_file, name) && !overwrite) {
    cli_abort(
      c(
        "Element {.val {name}} already exists in {.path {hdf5_file$path}}",
        "i" = "Set {.code overwrite = TRUE} to replace it"
      ),
      call = call
    )
  }

  parent <- dirname(name)
  if (parent != ".") {
    hdf5_ensure_group_path(hdf5_file, parent)
  }

  write_h5ad_element(value, hdf5_file, name, stop_on_error = TRUE, ...)

  invisible(NULL)
}

#' Write an element to a Zarr store
#'
#' @param store Path to an existing Zarr store
#' @param name Path of the element inside the store
#' @param value The value to write
#' @param overwrite Whether to replace an existing element
#' @param zarr_format The Zarr format of `store`
#' @param ... Passed to `write_zarr_element()`
#' @param call The calling environment, used in error messages
#'
#' @return `NULL` (invisibly)
#'
#' @noRd
.write_elem_zarr <- function(
  store,
  name,
  value,
  overwrite = FALSE,
  zarr_format,
  ...,
  call = rlang::caller_env()
) {
  name <- gsub("^/+|/+$", "", name)

  # TODO: Allow writing AnnData objects to the root
  if (name == "") {
    cli_abort("Writing to the root group is not supported yet", call = call)
  }

  if (zarr_path_exists(store, name) && !overwrite) {
    cli_abort(
      c(
        "Element {.val {name}} already exists in {.path {store}}",
        "i" = "Set {.code overwrite = TRUE} to replace it"
      ),
      call = call
    )
  }

  parent <- dirname(name)
  if (parent != "." && !zarr_path_exists(store, parent)) {
    Rarr::write_zarr_group(store, parent, zarr_version = zarr_format)
  }

  write_zarr_element(
    value,
    store,
    name,
    zarr_format = zarr_format,
    stop_on_error = TRUE,
    ...
  )

  # TODO: Update consolidated metadata once #492 is merged

  invisible(NULL)
}

#' Check element writing options
#'
#' @param allowed Character vector of supported argument names
#' @param ... The arguments to check
#' @param call The calling environment, used in error messages
#'
#' @return `NULL` (invisibly)
#'
#' @noRd
check_elem_dots <- function(allowed, ..., call = rlang::caller_env()) {
  dot_names <- names(rlang::list2(...)) %||% rep("", ...length())

  if (!all(nzchar(dot_names))) {
    cli_abort("All arguments in {.arg ...} must be named", call = call)
  }

  unknown <- setdiff(dot_names, allowed)
  if (length(unknown) > 0) {
    cli_abort(
      c(
        "Unsupported argument{?s} in {.arg ...}: {.arg {unknown}}",
        "i" = "Supported arguments are: {.arg {allowed}}"
      ),
      call = call
    )
  }

  invisible(NULL)
}
