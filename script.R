#'   library(SeuratObject)
#'   
#'   counts <- matrix(1:15, 3L, 5L)
#'   dimnames(counts) <- list(
#'     letters[1:3],
#'     LETTERS[1:5]
#'   )
#'   gene.metadata <- data.frame(
#'     row.names = LETTERS[1:5],
#'     gene = 1:5
#'   )
#'   obj <- CreateSeuratObject(counts, meta.data = gene.metadata)
#'   cell.metadata <- data.frame(
#'     row.names = letters[1:3],
#'     cell = 1:3
#'   )
#'   obj <- AddMetaData(obj, cell.metadata)
#'   