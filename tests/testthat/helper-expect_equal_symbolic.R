

expect_equal_symbolic <- function(object, expected){

  symbolic_difference <- Ryacas::yac_str(paste0("Simplify((", object,")",
                                                " - ",
                                                "(",expected,"))"))

  testthat::expect_equal(symbolic_difference, "0")
}
