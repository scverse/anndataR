#' S3 Methods for AbstractAnnData Objects
#'
#' These S3 methods provide standard R interfaces for AbstractAnnData objects,
#' making them behave like native R objects with familiar syntax.
#'
#' @name AbstractAnnData-s3methods
#' @param x An AbstractAnnData object
#' @param i Row indices (observations). Can be numeric, logical, or character.
#' @param j Column indices (variables). Can be numeric, logical, or character.
#' @param drop Ignored (for compatibility with generic)
#' @param ... Additional arguments passed to methods
#'
#' @details
#' The following S3 methods are available:
#' \itemize{
#'   \item \code{dim(x)}: Get dimensions (n_obs, n_vars)
#'   \item \code{nrow(x)}: Get number of observations
#'   \item \code{ncol(x)}: Get number of variables
#'   \item \code{rownames(x)}: Get observation names
#'   \item \code{colnames(x)}: Get variable names
#'   \item \code{x[i, j]}: Subset observations and/or variables
#' }
#'
#' @return
#' \itemize{
#'   \item \code{dim}: Numeric vector of length 2 (n_obs, n_vars)
#'   \item \code{nrow}, \code{ncol}: Integer count
#'   \item \code{rownames}, \code{colnames}: Character vector
#'   \item \code{[}: A AnnDataView object with the specified subset
#' }
#'
#' @examples
#' \dontrun{
#' # Create example data
#' ad <- generate_dataset(n_obs = 100, n_vars = 50)
#'
#' # Standard R methods work
#' dim(ad)
#' nrow(ad)
#' ncol(ad)
#' rownames(ad)
#' colnames(ad)
#'
#' # Subsetting creates AnnDataView
#' subset_ad <- ad[1:10, 1:5]
#' subset_ad <- ad[c(TRUE, FALSE), ]  # logical subsetting
#' subset_ad <- ad[c("obs_1", "obs_2"), c("var_1", "var_2")]  # name subsetting
#' }
NULL

#' @rdname AbstractAnnData-s3methods
#' @method dim AbstractAnnData
#' @export
dim.AbstractAnnData <- function(x) {
  x$shape()
}

#' @rdname AbstractAnnData-s3methods
#' @method nrow AbstractAnnData
#' @export
nrow.AbstractAnnData <- function(x) {
  x$n_obs()
}

#' @rdname AbstractAnnData-s3methods
#' @method ncol AbstractAnnData
#' @export
ncol.AbstractAnnData <- function(x) {
  x$n_vars()
}

#' @rdname AbstractAnnData-s3methods
#' @method rownames AbstractAnnData
#' @export
rownames.AbstractAnnData <- function(x) {
  x$obs_names
}

#' @rdname AbstractAnnData-s3methods
#' @method colnames AbstractAnnData
#' @export
colnames.AbstractAnnData <- function(x) {
  x$var_names
}

#' @rdname AbstractAnnData-s3methods
#' @method [ AbstractAnnData
#' @export
`[.AbstractAnnData` <- function(x, i, j, drop = TRUE, ...) {
  view <- AnnDataView$new(x)

  # Handle observation (row) subsetting
  if (!missing(i)) {
    if (is.logical(i)) {
      if (length(i) != x$n_obs()) {
        cli_abort(
          "Logical vector for observations must have length {x$n_obs()}"
        )
      }
    } else if (is.numeric(i)) {
      if (any(i < 1 | i > x$n_obs())) {
        cli_abort("Observation indices must be between 1 and {x$n_obs()}")
      }
      # Convert to logical
      logical_i <- rep(FALSE, x$n_obs())
      logical_i[i] <- TRUE
      i <- logical_i
    } else if (is.character(i)) {
      obs_names <- x$obs_names
      if (any(!i %in% obs_names)) {
        missing_names <- i[!i %in% obs_names]
        cli_abort(
          "Observation names not found: {paste(missing_names, collapse = ', ')}"
        )
      }
      # Convert to logical
      logical_i <- obs_names %in% i
      i <- logical_i
    }

    view$subset_obs(i)
  }

  # Handle variable (column) subsetting
  if (!missing(j)) {
    if (is.logical(j)) {
      if (length(j) != x$n_vars()) {
        cli_abort("Logical vector for variables must have length {x$n_vars()}")
      }
    } else if (is.numeric(j)) {
      if (any(j < 1 | j > x$n_vars())) {
        cli_abort("Variable indices must be between 1 and {x$n_vars()}")
      }
      # Convert to logical
      logical_j <- rep(FALSE, x$n_vars())
      logical_j[j] <- TRUE
      j <- logical_j
    } else if (is.character(j)) {
      var_names <- x$var_names
      if (any(!j %in% var_names)) {
        missing_names <- j[!j %in% var_names]
        cli_abort(
          "Variable names not found: {paste(missing_names, collapse = ', ')}"
        )
      }
      # Convert to logical
      logical_j <- var_names %in% j
      j <- logical_j
    }

    view$subset_var(j)
  }

  view
}
