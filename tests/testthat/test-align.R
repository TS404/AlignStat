context("Calculating identity score")

test_that("identity_score() works", {
  data(ref2)
  data(aln2)
  data(results)

  scores <- rcpp_align(ref2,aln2)

  expect_equal(class(scores),"matrix")
  expect_equal(scores[1,],results[1,])
  
})

