#' Check required packages
#'
#' Check that required packages are available and give a nice error message with
#' install instructions if not
#'
#' @param what A message stating what the packages are required for. Used at the
#'   start of the error message e.g. "{what} requires...".
#' @param requires Character vector of required package names
#' @param where Where to install the packages from. Either "CRAN", "Bioc", or "Python"
#' @param version Optional character vector with the minimum required version
#'   of each package in `requires`, `NA` for no minimum. Only used for R packages.
#'
#' @return `'TRUE` invisibly if all packages are available, otherwise calls
#'   [cli::cli_abort()]
#'
#' @importFrom rlang caller_env
#' @noRd
check_requires <- function(
  what,
  requires,
  where = c("CRAN", "Bioc", "Python"),
  version = NULL
) {
  where <- match.arg(where)

  if (where == "Python") {
    check_python_packages(what, requires)
  } else {
    check_r_packages(what, requires, where, version)
  }

  invisible(TRUE)
}

#' Check required Python packages
#'
#' @param what A message stating what the packages are required for
#' @param requires Character vector of required Python package names
#'
#' @importFrom rlang caller_env
#' @noRd
check_python_packages <- function(what, requires) {
  check_requires(what, "reticulate")

  is_available <- map_lgl(requires, reticulate::py_module_available)

  if (!all(is_available)) {
    missing <- requires[!is_available]
    # nolint start object_usage_linter
    missing_str <- format_package_list(missing)
    # nolint end object_usage_linter

    cli_abort(
      c(
        "{what} requires the Python {.pkg {missing}} package{?s}",
        "i" = paste(
          "To continue, require {cli::qty(missing)}{?it/them} using",
          "{.code reticulate::py_require({missing_str})} or install",
          "{cli::qty(missing)}{?it/them} them in your Python environment",
          "using another method"
        )
      ),
      call = caller_env()
    )
  }
}

#' Check required R packages
#'
#' @param what A message stating what the packages are required for
#' @param requires Character vector of required R package names
#' @param where Where to install the packages from. Either "CRAN" or "Bioc"
#' @param version Optional character vector with the minimum required version
#'   of each package in `requires`, `NA` for no minimum
#'
#' @importFrom rlang caller_env
#' @noRd
check_r_packages <- function(what, requires, where, version = NULL) {
  if (is.null(version)) {
    version <- rep(NA_character_, length(requires))
  }

  is_available <- map_lgl(seq_along(requires), function(i) {
    min_version <- if (is.na(version[[i]])) NULL else version[[i]]
    rlang::is_installed(requires[[i]], version = min_version)
  })

  if (!all(is_available)) {
    missing <- requires[!is_available]
    missing_version <- version[!is_available]
    # nolint start object_usage_linter
    missing_str <- format_package_list(missing)
    install_fun <- get_r_install_function(where)
    has_version <- !is.na(missing_version)
    version_str <- paste(
      missing[has_version],
      ">=",
      missing_version[has_version]
    )
    # nolint end object_usage_linter

    msg <- c(
      "{what} requires the {.pkg {missing}} package{?s}",
      "i" = paste(
        "To continue, install {cli::qty(missing)}{?it/them} using",
        "{.run {install_fun}({missing_str})}"
      )
    )
    if (any(has_version)) {
      msg <- c(
        msg,
        "i" = "{cli::qty(version_str)}Minimum version{?s} required: {.val {version_str}}"
      )
    }

    cli_abort(msg, call = caller_env())
  }
}

#' Get the version of a Python package
#'
#' @param package Name of the Python package
#'
#' @return The version as a [package_version()] object. Pre-release and local
#'   suffixes (e.g. `0.13.0rc1`) are dropped.
#' @noRd
py_package_version <- function(package) {
  metadata <- reticulate::import("importlib.metadata", convert = TRUE)
  version <- metadata$version(package)
  package_version(sub("^([0-9]+(\\.[0-9]+)*).*$", "\\1", version))
}

#' Format package list for installation command
#'
#' @param packages Character vector of package names
#' @return Formatted string for use in installation functions
#' @noRd
format_package_list <- function(packages) {
  missing_str <- paste0("\"", paste(packages, collapse = "\", \""), "\"")
  if (length(packages) > 1) {
    missing_str <- paste0("c(", missing_str, ")")
  }
  missing_str
}

#' Get R package installation function
#'
#' @param where Where to install from. Either "CRAN" or "Bioc"
#' @return Installation function name as string
#' @noRd
get_r_install_function <- function(where) {
  # nolint start object_usage_linter
  switch(
    where,
    CRAN = "install.packages",
    Bioc = "install.packages(\"BiocManager\"); BiocManager::install"
  )
  # nolint end object_usage_linter
}
