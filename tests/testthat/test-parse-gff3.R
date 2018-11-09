context("Read GFF3 files.")
library(codonfriend)
library(rtracklayer)

library(GenomicRanges)


# This is our dummy file that we will use.
tmp_gff1 <- tempfile(fileext = ".gff3")
gff1_text <- "##gff-version 3
seq1\tannot\tgene\t1\t5\t0\t+\t.\tID=id1;Name=name1
seq1\tannot\tmRNA\t1\t5\t0\t+\t.\tID=id1.t1;Parent=id1
"

setup({
  # Write the dummy gff to a tempfile.
  writeLines(text = gff1_text, con = tmp_gff1)
})

teardown({
  # Deletes the tempfiles
  unlink(tmp_gff1)
})

test_that("read_gff parses file", {
  gff <- read_gff(path = tmp_gff1)

  # Two lines in GFF file
  expect_equal(length(gff), 2)

  # Starts are correct
  expect_equal(start(ranges(gff)), c(1, 1))

  # Ends are correct
  expect_equal(start(ranges(gff)), c(5, 5))

  # Specified ids are correct.
  expect_equal(mcols(gff)[, "ID"], c("id1", "id1.t1"))
  
  # Feature types are correct.
  expect_equal(mcols(gff)[, "type"], c("gene", "mRNA"))
})