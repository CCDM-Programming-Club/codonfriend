context("Join CAI with GRanges object.")
library(codonfriend)
library(GenomicRanges)


test_that("join_cai can join simple files.", {

  gr1 <- GRanges(
    seqnames = rep("C1", 4),
    ranges = IRanges(
      start = c(1, 11, 21, 31),
      end = c(6, 16, 26, 36),
      names = c("A", "B", "C", "D")
    ),
    strand = rep("+", 4),
    ID=c("g1", "g2", "g3", "g4")
  )
  
  cai1 <- data.frame(
    seqid = c("g1", "g2", "g3", "g4"),
    cai = c(1.0, 0.1, 0.5, 0.3)
  )
  
  result <- join_cai(gr1, cai1)

  # Did we keep all of the rows and add the new column?
  expect_equal(dim(mcols(result)), c(4, 2))

  # Are the columns the names and order that we'd expect?
  expect_equal(
    names(mcols(result)),
    c("ID", "cai")
  )

  # Do the CAI values get added?
  expect_equal(
    mcols(result[c("B", "D", "C", "A")])[["cai"]],
    c(0.1, 0.3, 0.5, 1.0)
  )
})


test_that("join_cai can handle extra columns and odd ID values.", {

  gr2 <- GRanges(
    seqnames = rep("C1", 7),
    ranges = IRanges(
      start = c(1, 2, 2, 11, 12, 21, 22),
      end = c(6, 6, 6, 16, 16, 26, 26),
      names = c("A", "B", "C", "D", "E", "F", "G")
    ),
    strand = rep("+", 7),
    ID = c("g1", "g1.t1", "g1.t1.e1", "g2", "g2.t1", "g3", NA),
    type = c("gene", "mRNA", "CDS", "gene", "mRNA", "gene", "CDS")
  )
  
  cai2 <- data.frame(
    seqid = c("g1.t1", "g2.t1", "g3"),
    cai = c(1.0, 0.1, 0.5),
    class = c("H", "I", "H")
  )
  
  result <- join_cai(gr2, cai2)

  # Did we keep all of the rows and add the new columns?
  expect_equal(dim(mcols(result)), c(7, 4))
  
  # Are the columns the names and order that we'd expect?
  expect_equal(
    names(mcols(result)),
    c("ID", "type", "cai", "class")
  )

  # Do the CAI values get added, with NAs filling other values?
  expect_equal(
    mcols(result[c("B", "F", "E", "A", "C", "D", "G")])[["cai"]],
    c(1.0, 0.5, 0.1, NA, NA, NA, NA)
  )

  # Do the class values get added?
  expect_equal(
    mcols(result[c("B", "F", "E", "A", "C", "D", "G")])[["class"]],
    as.factor(c("H", "H", "I", NA, NA, NA, NA))
  )
})
