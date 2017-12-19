context("basic functionality")
test_that("we can do something", {

  x <- saucy::saucy(system.file("extdata", "graphfile", package="saucy"))
  expect_is(x, "saucy")

  y <- saucy::saucy(system.file("extdata", "graphfile2", package="saucy"))
  expect_is(y, "saucy")

  z <- saucy::shatter(system.file("extdata", "example.cnf", package="saucy"))
  expect_is(z, "shatter")

})



