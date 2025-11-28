test_that("get_generator_types() works", {
  types <- get_generator_types()

  expect_type(types, "list")
  expect_true(all(.anndata_slots %in% names(types)))
  expect_true(all(sapply(types, is.character)))
})

test_that("get_generator_types(example = TRUE) works", {
  types <- get_generator_types(example = TRUE)

  expect_type(types, "list")
  expect_true(all(.anndata_slots %in% names(types)))
  expect_true(all(sapply(types, is.character)))

  for (slot in names(types)) {
    expect_true(all(types[[slot]] %in% .example_generator_types[[slot]]))
  }
})

test_that("get_generator_types(slot = ...) works", {
  for (slot in .anndata_slots) {
    types <- get_generator_types(slot = slot)

    expect_type(types, "character")
    expect_true(all(types %in% .generator_types[[slot]]))
  }
})
