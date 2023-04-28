test_that("AnnData_to_SingleCellExperiment() works", {
    ad <- InMemoryAnnData$new(
        X = matrix(1:5, 3L, 5L),
        obs = data.frame(cell = 1:3, row.names = LETTERS[1:3]),
        var = data.frame(gene = 1:5, row.names = letters[1:5])
    )

    expect_no_error(sce <- to_SingleCellExperiment(ad))
    expect_true(validObject(sce))
    expect_identical(dim(sce), rev(ad$shape()))

    ad_dimnames <- list(rownames(ad$var), rownames(ad$obs))
    expect_identical(dimnames(sce), ad_dimnames)
    expect_identical(as.data.frame(SummarizedExperiment::rowData(sce)), ad$var)
    expect_identical(as.data.frame(SummarizedExperiment::colData(sce)), ad$obs)

    expect_identical(
        SummarizedExperiment::assay(sce, withDimnames = FALSE),
        t(ad$X)
    )
})
