#' Calculate coverage log-ratio of tumor vs. normal
#' 
#' This function is automatically called by \code{\link{runAbsoluteCN}} when
#' normal and tumor coverage are provided (and not a segmentation file or
#' target-level log-ratios). This function is therefore normally not called by
#' the user.
#' 
#' 
#' @param normal Normal coverage read in by the \code{\link{readCoverageFile}}
#' function.
#' @param tumor Tumor coverage read in by the \code{\link{readCoverageFile}}
#' function.
#' @return \code{numeric(length(tumor))}, tumor vs. normal copy number log-ratios
#' for all targets.
#' @author Markus Riester
#' @examples
#' 
#' normal.coverage.file <- system.file("extdata", "example_normal.txt.gz", 
#'     package = "PureCN")
#' tumor.coverage.file <- system.file("extdata", "example_tumor.txt.gz", 
#'     package = "PureCN")
#' normal <- readCoverageFile(normal.coverage.file)
#' tumor <- readCoverageFile(tumor.coverage.file)
#' log.ratio <- calculateLogRatio(normal, tumor)
#' 
#' @export calculateLogRatio
calculateLogRatio <- function(normal, tumor) {
    # make sure that normal and tumor align
    if (!identical(as.character(normal), as.character(tumor))) {
        .stopUserError("Interval files in normal and tumor different.")
    }
    if (is.null(tumor$on.target)) tumor$on.target <- TRUE

    avgCovTumor <- mean(tumor$average.coverage[tumor$on.target], na.rm=TRUE)
    avgCovNormal <- mean(normal$average.coverage[tumor$on.target], na.rm=TRUE)

    flog.info("Mean target coverages: %.0fX (tumor) %.0fX (normal).", 
        avgCovTumor, avgCovNormal)
    if (avgCovNormal/avgCovTumor < 0.25 || avgCovNormal/avgCovTumor > 4) {
        flog.warn("Large difference in coverage of tumor and normal.")
    }    
    
    tumor$log.ratio <- 0.

    for (on.target in c(FALSE, TRUE)) {
        idx <- tumor$on.target==on.target
        if (!sum(idx)) next
        total.cov.normal <- sum(as.numeric(normal[idx]$coverage), na.rm = TRUE)
        total.cov.tumor <- sum(as.numeric(tumor[idx]$coverage), na.rm = TRUE)

        log.ratio <- log2(tumor[idx]$average.coverage/normal[idx]$average.coverage) + 
                     log2(total.cov.normal/total.cov.tumor)
        tumor[idx]$log.ratio <- .calibrate_log_ratio(log.ratio, tumor[idx])
    }
    if (!all(tumor$on.target)) {
        # try to align the off-target and on-target log-ratios better
        tumor$log.ratio <- .calibrate_off_target_log_ratio(tumor)
    }    
    tumor$log.ratio
}

.calibrate_log_ratio <- function(log.ratio, granges) {
    idxFinite <- is.finite(log.ratio)
    if (!sum(idxFinite)) {
        .stopUserError("No finite intervals.")
    }
    mean.log.ratio <- weighted.mean(log.ratio[idxFinite], 
        width(granges)[idxFinite])
# calibrate
    flog.debug("Calibrating %i log-ratios by %f.",
               sum(idxFinite), mean.log.ratio)
    return(log.ratio - mean.log.ratio)
}

.calibrate_off_target_log_ratio <- function(granges) {
    idx <- granges$on.target
    g1 <- granges[idx]
    g2 <- granges[!idx]
    nr <- nearest(g1,g2)
    d2 <- median(g1$log.ratio - g2$log.ratio[nr], na.rm = TRUE) / 2
    if (d2 > 0.1) {
        flog.warn("Large potential mis-calibration of on- and off-target log2 ratios: %.2f", d2)
    }
    granges$log.ratio[idx] <- granges$log.ratio[idx] - d2
    granges$log.ratio[!idx] <- granges$log.ratio[!idx] + d2
    return(granges$log.ratio)
}
