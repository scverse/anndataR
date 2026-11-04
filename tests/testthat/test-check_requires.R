test_that("check_requires passes for an installed package", {
  expect_true(check_requires("Test", "rlang"))
})

test_that("check_requires errors for a missing package", {
  cnd <- expect_error(check_requires("Test", "notapackage"))
  expect_match(conditionMessage(cnd), "notapackage")
  expect_no_match(conditionMessage(cnd), "Minimum version")
})

test_that("check_requires passes when the installed version is new enough", {
  expect_true(check_requires("Test", "rlang", version = "0.0.1"))
})

test_that("check_requires errors when the installed version is too old", {
  cnd <- expect_error(check_requires("Test", "rlang", version = "999.0.0"))
  expect_match(conditionMessage(cnd), "rlang")
  expect_match(conditionMessage(cnd), "999.0.0")
})
