test_that("basic one person likelihood works", {

  require(pedtools)

  het <- singleton() |> addMarker(geno="a/b")
  expect_equal_symbolic(sLikelihood(het), "a*b+a*b")

  hom <- singleton() |> addMarker(geno="a/a")
  expect_equal_symbolic(sLikelihood(hom), "a*a")

})
