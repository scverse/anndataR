#' Check required packages
#'
#' Check that required packages are available and give a nice error message with
#' install instructions if not
#'
#' @param what A message stating what the packages are required for. Used at the
#'   start of the error message e.g. "{what} requires...".
#' @param requires Character vector of required package names
#' @param where Where to install the packages from. Either "CRAN" or "Bioc"
#'
#' @return `'TRUE` invisibly if all packages are available, otherwise calls
#'   [cli::cli_abort()]
#'
#' @importFrom rlang caller_env
#' @noRd
check_requires <- function(what, requires, where = c("CRAN", "Bioc")) {
  where <- match.arg(where)

  is_available <- map_lgl(requires, requireNamespace, quietly = TRUE)

  if (any(!is_available)) {
    missing <- requires[!is_available]

    # nolint start object_usage_linter
    missing_str <- paste0("\"", paste(missing, collapse = "\", \""), "\"")
    if (length(missing) > 1) {
      missing_str <- paste0("c(", missing_str, ")")
    }
    fun <- switch(
      where,
      CRAN = "install.packages",
      Bioc = "install.packages(\"BiocManager\"); BiocManager::install"
    )
    # nolint end object_usage_linter

    cli_abort(
      c(
        "{what} requires the {.pkg {missing}} package{?s}",
        "i" = paste(
          "To continue, install {cli::qty(missing)}{?it/them} using",
          "{.code {fun}({missing_str})}"
        )
      ),
      call = caller_env()
    )
  }

  invisible(TRUE)
}
