
#' Parse a file from EMBOSS cai program.
#'
#' @param path A string containing the path to a cai file.
#' @return A dataframe containing the parsed data with column headings (seqid, cai).
#' @examples
#' # Load some sample data included with the package.
#' cai_path <- system.file("extdata", "lmac.cai", package = "codonfriend")
#' # Parse that file.
#' cai_table <- read_cai(cai_path)
#' head(cai_table)
#' 
#' @importFrom stringr str_trim
#' @importFrom stringr str_match
#' @export
read_cai <- function(path) {
  lines <- readLines(path)

  # Filters out empty lines.
  # Could also handle comment lines here
  clean_lines <- subset(
    lines,
    (str_trim(lines) != "")
  )

  # Extract the important stuff using a regular expression
  reg <- "Sequence: (.+) CAI: (.+)"
  matches <- str_match(clean_lines, reg)

  data <- data.frame(
     seqid = str_trim(matches[, 2]),
     cai = as.numeric(matches[, 3]),
     stringsAsFactors = FALSE
  )
  return(data)
}


#' Read many files and combine them.
#' 
#' @param paths A vector of files to parse.
#' @param FUN The function to use to parse the files. Default `read_cai`.
#' @param colname What to call the column added.
#' @return A concatenated dataframe with filenames added as a new column.
#' @examples
#' 
#' paths <- list.files("dir", full.names = TRUE)
#' # or
#' paths <- Sys.glob("dir/*.cai")
#' df <- read_many(paths, FUN=read_cai, colname="file")
#' head(df)
#' 
#' @export
read_many <- function(paths, FUN=read_cai, colname="file") {
  do.call(
    rbind,
    args = lapply(
      paths,
      FUN = function(f) {
        t <- FUN(f)
        t[colname] <- basename(f)
        return(t)
      }
    )
  )
}


#' Read a GFF3 file as a GRanges object..
#' 
#' @param path The path to the GFF file.
#' @return A GenomicRanges `GRanges` object.
#' @examples
#' 
#' path <- system.file("extdata", "Lepmu1.gff3", package = "codonfriend")
#' gff <- read_gff3(path)
#' head(gff)
#' 
#' @importFrom rtracklayer
#' @importFrom GenomicRanges GRanges
#' @export
read_gff3 <- function(path) {
  # Add your code here! Paula
  gff.file <- read.gff(path)
  # Import files by adding to docstring like in `parse_cai`. E.g.
  # #' @importFrom GenomicRanges GRanges
  gff3 <- readGFF(path)
  return(gff3)

}
