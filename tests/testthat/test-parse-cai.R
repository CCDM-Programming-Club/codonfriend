context("Parse cai successfully.")
library(codonfriend)

test_that("read_cai Parses well formatted file", {
  string <- "Sequence: gene1 CAI: 0.1
Sequence: gene2 CAI: 0.2"

  # NB textConnection creates a file-like object for readLines to work on.
  con <- textConnection(string)
  table <- read_cai(con)

  expect_equal(dim(table), c(2, 2))
  expect_equal(names(table), c("seqid", "cai"))
  expect_equal(table[, "seqid"], c("gene1", "gene2"))
  expect_equal(table[, "cai"], c(0.1, 0.2))
  expect_equal(sapply(table, class), c(seqid = "character", cai = "numeric"))
})

test_that("read_cai Strips trailing empty lines", {
  # Here the trailing newline could add an empty row.
  string <- "Sequence: gene1 CAI: 0.1
Sequence: gene2 CAI: 0.2
"

  con <- textConnection(string)
  table <- read_cai(con)

  expect_equal(dim(table), c(2, 2))
})

test_that("read_cai Handles odd seqids", {
  string <- "Sequence: gene1#garbage CAI: 0.1
Sequence: gene2 with spaces CAI: 0.2
Sequence: RS|001_05-gene.oo@ CAI: 0.3
Sequence: *& ^%$#@!/? CAI: 0.4
Sequence: a      CAI: 0.4"

  con <- textConnection(string)
  table <- read_cai(con)

  expect_equal(dim(table), c(5, 2))
  expect_equal(
    table[, "seqid"],
    c("gene1#garbage", "gene2 with spaces",
      "RS|001_05-gene.oo@", "*& ^%$#@!/?", "a")
  )
})

test_that("read_cai Invalid CAIs as NA", {
  string <- "Sequence: gene1 CAI: 0.1
Sequence: gene2 CAI: .
Sequence: gene3 CAI: NA"

  con <- textConnection(string)
  expect_warning(read_cai(con), "NAs introduced by coercion")

  # NB con is consumed by read_cai, need second object
  con <- textConnection(string)
  table <- suppressWarnings(read_cai(con))

  expect_equal(dim(table), c(3, 2))
  expect_equal(table[, "cai"], c(0.1, NA, NA))
})

test_that("read_cai Empty file as zero rows", {
  string <- ""

  con <- textConnection(string)
  table <- read_cai(con)

  expect_equal(dim(table), c(0, 2))
})
