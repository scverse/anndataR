# helper function to skip tests if h5diff is not available
skip_if_no_h5diff <- function() {
  testthat::skip_if(
    !system("which h5diff", ignore.stdout = TRUE),
    message = "h5diff not available for testing"
  )
}
