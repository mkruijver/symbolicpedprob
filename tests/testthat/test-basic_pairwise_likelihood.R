test_that("basic pairwise likelihood works", {

  require(pedtools)
  x <- nuclearPed(nch = 2) |>
    addMarker(geno = c(NA, NA, "a/b", "c/d"))

  lik <- sLikelihood(x)
  expected_lik <- "a*b*c*d"

  expect_equal_symbolic(lik, expected_lik)
})
