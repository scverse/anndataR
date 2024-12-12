# helper function to skip tests if h5diff is not available
skip_if_no_h5diff <- function() {
  testthat::skip_if_not({
    s <- system2(command = "which", args = "h5diff", stdout = TRUE, stderr = TRUE)
    is.null(attr(s, "status"))
  }, message = "h5diff not available for testing")
}
