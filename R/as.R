#' Convert an `AnnData` to an `InMemoryAnnData`
#'
#' Convert another `AnnData` object to an [`InMemoryAnnData`] object
#'
#' @param adata An `AnnData` object to be converted to [`InMemoryAnnData`]
#'
#' @return An [`InMemoryAnnData`] object with the same data as the input
#'   `AnnData` object
#' @name as_InMemoryAnnData
#'
#' @family object converters
#'
#' @examples
#' ad <- AnnData(
#'   X = matrix(1:5, 3L, 5L),
#'   layers = list(
#'     A = matrix(5:1, 3L, 5L),
#'     B = matrix(letters[1:5], 3L, 5L)
#'   ),
#'   obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
#'   var = data.frame(row.names = letters[1:5], gene = 1:5)
#' )
#' ad$as_InMemoryAnnData()
NULL
