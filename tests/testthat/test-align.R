context("Calculating identity score")

test_that("rpp_align() produces correct results", {
  data(ref2)
  data(com2)
  data(results)
  data(means)
  data(cat2)

  scores <- rcpp_align(ref2,com2)

  expect_equal(class(scores),"list")
  sresults = scores$results
  expect_equal(sresults[1,],results[1,])
  smeans = scores$means
  expect_equal(smeans,means)
  
  expect_equal(scores$cat,cat2)
})


test_that("align_alignments() produces correct results", {
  data(ref)
  data(com)
  data(results)
  data(means)
  
  scores <- align_alignments(ref,com)
  
  expect_equal(class(scores),"list")
  sresults = scores$results
  expect_equal(sresults,results)
  smeans = scores$means
  expect_equal(smeans,means)
})


test_that("prepare_alignment_matrix() produces correct outputs",{
  data(ref)
  data(ref2)
  r2 <- prepare_alignment_matrix(ref)
  
  expect_equal(ref2,r2)
  
})

context("Performance")

test_that("align_alignments() is fast enough",{
  data(ref)
  data(com)
  timing <- system.time(align_alignments(ref,com))
  
  cat("\nAlign alignments took: ",timing)
  expect_less_than(timing[1],1)
})


