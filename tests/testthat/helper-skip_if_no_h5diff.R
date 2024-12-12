# helper function to skip tests if h5diff is not available
skip_if_no_h5diff <- function() {
  testthat::skip_if(
    tryCatch({
      system2("which h5diff")
      FALSE
    }, error = function(e) TRUE
    ),
    message = "h5diff not available for testing"
  )
}
