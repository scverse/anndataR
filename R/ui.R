#' Style a vector
#'
#' Style a vector for use in **{cli}** calls
#'
#' @param x The vector to style
#' @param last The separator before the last element in the vector
#' @param trunc The maxium number of elements to show
#' @param trunc_style The truncation style to use, either "both-ends"
#'   (1, 2, ..., n-1, 1) or "head" (1, 2, 3, 4, ...)
#' @param wrap Whether to wrap the output to appear as a vector, if `TRUE` then
#' "c(1, 2, 3)" else "1, 2, 3"
#'
#' @returns A styled string that can be included directly using "{syle_vec(...)}"
#' @noRd
style_vec <- function(
  x,
  last = ", ",
  trunc = 20L,
  trunc_style = c("both-ends", "head"),
  wrap = FALSE
) {
  trunc_style <- match.arg(trunc_style)

  style <- list(
    "vec-last" = last,
    "vec-trunc" = trunc,
    "vec-trunc-style" = trunc_style
  )

  x <- cli::cli_vec(x, style)

  if (isTRUE(wrap)) {
    x_str <- "c({.val {x}})"
  } else {
    x_str <- "{.val {x}}"
  }

  cli::cli_fmt(cli::cli_text(x_str))
}
