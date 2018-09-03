library("stringr")

read.cai <- function(path) {
  lines <- readLines(path)
  reg <- "Sequence: (.+) CAI: (.+)"
  matches <- str_match(lines, reg)
  
  data <- data.frame(
     seqid = as.character(matches[,2]),
     CAI = as.numeric(matches[,3]),
     stringsAsFactors = FALSE
  )
  return(data)
}
