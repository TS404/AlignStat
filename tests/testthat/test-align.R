context("Calculating identity score")

test_that("rpp_align() produces correct results", {
  data(ref2)
  data(aln2)
  data(results)
  data(means)

  scores <- rcpp_align(ref2,aln2)

  expect_equal(class(scores),"list")
  sresults = scores$results
  expect_equal(sresults[1,],results[1,])
  smeans = scores$means
  expect_equal(smeans,means)
})


test_that("align_alignments() produces correct results", {
  data(ref)
  data(aln)
  data(results)
  data(means)
  
  scores <- align_alignments(ref,aln)
  
  expect_equal(class(scores),"list")
  sresults = scores$results
  expect_equal(sresults,results)
  smeans = scores$means
  expect_equal(smeans,means)
})

