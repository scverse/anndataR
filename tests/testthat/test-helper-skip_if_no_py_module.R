skip_if_not_installed("reticulate")

test_that("skip_if_no_py_module() skips when a module is unavailable", {
  # mock reticulate so that no Python session is started
  local_mocked_bindings(
    py_module_available = function(module) FALSE,
    import = function(module, ...) stop("No module named '", module, "'"),
    .package = "reticulate"
  )
  withr::local_envvar(ANNDATAR_TEST_REQUIRE_PYTHON = NA)

  expect_condition(
    skip_if_no_py_module("not_a_module", requirement = NULL),
    class = "skip"
  )
})

test_that("skip_if_no_py_module() errors when ANNDATAR_TEST_REQUIRE_PYTHON is set", {
  # mock reticulate so that no Python session is started
  local_mocked_bindings(
    py_module_available = function(module) FALSE,
    import = function(module, ...) stop("No module named '", module, "'"),
    .package = "reticulate"
  )
  withr::local_envvar(ANNDATAR_TEST_REQUIRE_PYTHON = "true")

  # catch the condition explicitly, a skip would otherwise escape expect_error()
  cnd <- rlang::catch_cnd(skip_if_no_py_module(
    "not_a_module",
    requirement = NULL
  ))
  expect_s3_class(cnd, "error")
  expect_false(inherits(cnd, "skip"))
  expect_match(conditionMessage(cnd), "No module named 'not_a_module'")
})
