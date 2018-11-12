#' Add CAI information to the metadata table of a GRanges object.
#'
#' @param gr A `GRanges` object corresponding to the genome the cais were calculated from.
#' @param cai A dataframe containing the columns, `seqid` and any number of other columns to be added e.g. `cai`.
#' @return An updated Granges object.
#' @examples
#' # Load some sample data included with the package.
#' cai_path <- system.file("extdata", "Lepmu1.cai", package = "codonfriend")
#' # Parse that file.
#' cai_table <- read_cai(cai_path)
#' 
#' gff_path <- system.file("extdata", "Lepmu1.gff3", package = "codonfriend")
#' gff <- rtracklayer::readGFFAsGRanges(gff_path)
#' 
#' gr <- join_cai(gff, cai)
#' head(gr)
#' 
#' @importFrom GenomicRanges mcols
#' @export
join_cai <- function(gr, cai) {
  meta <- mcols(gr, use.names = TRUE)
  meta$tmp_order <- 1:nrow(meta)
  merged <- merge(
    meta,
    cai,
    by.x="ID",
    by.y="seqid",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  merged <- merged[order(merged$tmp_order),]
  merged$tmp_order <- NULL
  mcols(gr) <- merged
  return(gr)
}