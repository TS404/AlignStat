context("Calculating identity score")

# test_that("rpp_align() produces correct results", {
#   prepared_ref <- readRDS("prepared_ref.rda")
#   prepared_com <- readRDS("prepared_com.rda")
#   results <- readRDS("results.rda")
#   means <- readRDS("means.rda")
# 
#   scores <- rcpp_align(prepared_ref,prepared_com)
# 
#   expect_equal(class(scores),"list")
#   sresults = scores$results
#   expect_equal(sresults[1,],results[1,])
#   smeans = scores$means
#   expect_equal(smeans,means)
# 
# })


test_that("compare_alignments() produces a list", {
  data(reference_alignment)
  data(comparison_alignment)

  scores <- compare_alignments(reference_alignment,comparison_alignment)

  expect_equal(class(scores),"list")
})

test_that("compare_alignments() handles non identical labels if sequences are identical", {
  data(reference_alignment)
  data(comparison_alignment)
  names(reference_alignment)[1]  <- "Wrong"
  
  expect_silent(compare_alignments(reference_alignment,comparison_alignment))
})

test_that("compare_alignments() errors if sequences are not identical", {
  data(reference_alignment)
  data(comparison_alignment)
  
  reference_alignment[[1]][1]="s"

  expect_error(compare_alignments(reference_alignment,comparison_alignment))
})

test_that("compare_alignments() works with file input", {
  file_a <- "AlignmentA.fasta"
  file_b <- "AlignmentB.fasta"
  scores <- compare_alignments(file_a,file_b)

  expect_equal(class(scores),"list")
})

test_that("compare_alignments() works with clustal and fasta input", {
  file_a <- "AlignmentA.fasta"
  file_b <- "AlignmentB.clustal"
  scores <- compare_alignments(file_a,file_b)
  
  expect_equal(class(scores),"list")
  expect_equal(scores$comlen,99)
})

test_that("compare_alignments() works with Grammicidins example", {
  file_a <- "GrammicidinsClustal.fasta"
  file_b <- "GrammicidinsMuscle.fasta"
  scores <- compare_alignments(file_a,file_b)
  
  expect_equal(class(scores),"list")
  expect_equal(scores$comlen,29)
})

