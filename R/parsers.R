
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
read_cai <- function(path) {
  lines <- readLines(path)

  # Filters out empty lines.
  # Could also handle comment lines here
  clean_lines <- subset(
    lines,
    (stringr::str_trim(lines) != "")
  )

  # Extract the important stuff using a regular expression
  reg <- "Sequence: (.+) CAI: (.+)"
  matches <- stringr::str_match(clean_lines, reg)

  data <- data.frame(
     seqid = stringr::str_trim(matches[,2]),
     cai = as.numeric(matches[,3]),
     stringsAsFactors = FALSE
  )
  return(data)
}
