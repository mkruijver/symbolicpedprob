
test_that("basic one person likelihood works for inbred person", {
  require(pedtools)

  x = cousinPed(1, child = T) |>
    addMarker(`9` = "a/a", alleles = c("a", "b"))

  # verify against known case
  x_expected <- "1/16*a + 15/16*a^2"
  expect_equal_symbolic(sLikelihood(x), x_expected)

  # verify against numerical calculation
  expect_equal(eval(Ryacas::yac_expr(sLikelihood(x)), list(a = 0.5)),
               pedprobr::likelihood(x))


  # another inbred pedigree
  y = cousinPed(2, child = TRUE)|>
    addMarker(`13` = "a/b", alleles = c("a", "b", "c", "d"))

  expect_equal(eval(Ryacas::yac_expr(sLikelihood(y)), list(a = 0.25, b = 0.25)),
               pedprobr::likelihood(y))
}
)
