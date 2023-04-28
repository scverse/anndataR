test_that("AnnData_to_SingleCellExperiment() works", {
    ad <- InMemoryAnnData$new(
        X = matrix(1:5, 3L, 5L),
        obs = data.frame(cell = 1:3),
        var = data.frame(gene = 1:5)
    )

    expect_no_error(sce <- to_SingleCellExperiment(ad))
    expect_true(validObject(sce))
    expect_identical(dim(sce), rev(ad$shape()))

    ad_dimnames <- list(rownames(ad$var), rownames(ad$obs))
    expect_identical(dimnames(sce), ad_dimnames)
    df <- as.data.frame(SummarizedExperiment::rowData(sce))
    rownames(df) <- NULL
    expect_identical(df, ad$var)
    df <- as.data.frame(SummarizedExperiment::colData(sce))
    rownames(df) <- NULL
    expect_identical(df, ad$obs)

    expect_identical(
        SummarizedExperiment::assay(sce, withDimnames = FALSE),
        t(ad$X)
    )
})
