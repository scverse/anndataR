test_that("example_data works", {
  
  classes <- c("InMemoryAnnData","HDF5AnnData","Seurat")
  for(x in classes){
    obj <- example_data(output_class = x)
    testthat::expect_true(methods::is(obj,x))
  } 
})
