#' Read interval file
#' 
#' Read file containing coordinates of on- and off-target intervals
#' generated by \code{\link{preprocessIntervals}}. 
#' 
#' @param interval.file A mapping file that assigns GC content and gene symbols
#' to each exon in the coverage files. Used for generating gene-level calls.
#' First column in format CHR:START-END. Second column GC content (0 to 1).
#' Third column gene symbol. This file is generated with the
#' \code{\link{preprocessIntervals}} function.
#' @param strict Error out with missing columns
#' @param verbose Verbose output
#' @return A \code{GRanges} object with the parsed intervals.
#' @author Markus Riester
#' @examples
#' 
#' interval.file <- system.file("extdata", "example_intervals.txt", 
#'     package = "PureCN")
#' x <- readIntervalFile(interval.file)
#' 
#' @export readIntervalFile
readIntervalFile <- function(interval.file, strict = TRUE, verbose = TRUE) {
    con <- file(interval.file, open = "r")
    header <- .parseGATKHeader(con)
    intervals <- read.delim(con, header = FALSE, stringsAsFactors = FALSE)
    colnames(intervals) <- strsplit(header$last_line, "\t")[[1]]
    close(con)
    if (is.null(intervals$gc_bias) && strict) {
        .stopUserError("No gc_bias column in interval.file.")
    }    
    if (is.null(intervals$Gene)) {
        if (verbose) flog.info("No Gene column in interval.file. You won't get gene-level calls.")
        intervals$Gene <- "."
    }
    if (is.null(intervals$on_target)) {
        if (verbose) flog.info("No on_target column in interval.file. Recreate this file with IntervalFile.R.")
        intervals$on_target <- TRUE
    }
    if (is.null(intervals$mappability)) {
        if (verbose) flog.info("No mappability column in interval.file.")
        intervals$mappability <- 1
    }
    if (is.null(intervals$reptiming)) {
        if (verbose) flog.info("No reptiming column in interval.file.")
        intervals$reptiming <- NA
    }
    
    gr <- GRanges(intervals[, 1], ranges = NULL, strand = NULL, intervals[, -1])
    gr <- sort(sortSeqlevels(gr))
    # TODO cleanup
    gr$on.target <- gr$on_target
    gr$on_target <- NULL

    if (length(header$sl)) {
        header$sl <- sapply(header$sl, as.numeric)
        seqlengths(gr) <- header$sl[names(seqlengths(gr))]
    }
    return(gr)
}
