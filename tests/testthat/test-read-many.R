context("Read multiple dataframe-like files with read_many successfully.")
library(codonfriend)

tmp_csv1 <- tempfile(fileext = ".csv")
tmp_csv2 <- tempfile(fileext = ".csv")
tmp_csv3 <- tempfile(fileext = ".csv")
tmp_csv4 <- tempfile(fileext = ".csv") # Empty

setup({
  writeLines(con = tmp_csv1, text = c("a,b,c,d", "e,f,g,h"))
  writeLines(con = tmp_csv2, text = c("i,j,k,l", "m,n,o,p"))
  writeLines(con = tmp_csv3, text = c("q,r,s,t", "u,v,w,x", "y,z,,"))
  writeLines(con = tmp_csv4, text = "")
})

teardown({
  unlink(tmp_csv1)
  unlink(tmp_csv2)
  unlink(tmp_csv3)
  unlink(tmp_csv4)
})



test_that("read_many parses many (3) files", {
  files <- c(tmp_csv1, tmp_csv2, tmp_csv3)
  table <- read_many(paths = files, FUN = function(x) read.csv(x, header = FALSE, stringsAsFactors = FALSE))
  
  # Dim is row-first
  expect_equal(dim(table), c(7, 5))
  
  expect("file" %in% names(table))
  expect_equal(table[, 1], c("a", "e", "i", "m", "q", "u", "y"))
})

test_that("read_many parses single file", {
  table <- read_many(paths = tmp_csv1, FUN = function(x) read.csv(x, header = FALSE, stringsAsFactors = FALSE))
  
  expect_equal(dim(table), c(2, 5))
  expect_equal(table[, 1], c("a", "e"))
})