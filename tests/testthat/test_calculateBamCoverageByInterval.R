context("calculateBamCoverageByInterval")

output.file <- tempfile(fileext = ".txt")

test_that("Coverage from test BAM file matches", {
    bam.file <- system.file("extdata", "ex1.bam", package = "PureCN", 
        mustWork = TRUE)
    interval.file <- system.file("extdata", "ex1_intervals.txt", 
        package = "PureCN", mustWork = TRUE)
    coverage <- calculateBamCoverageByInterval(bam.file = bam.file, 
        interval.file = interval.file, output.file = output.file)
    expect_equal(coverage$average.coverage, c(20.95205, 43.78357, 
        21.29271), tolerance = 0.01)
    expect_equal(coverage$counts, c(610, 1158, 636), tolerance = 0.01)
    expect_equal(unlist(coverage$duplication.rate), rep(0, 3), check.names = FALSE)
}) 

test_that("Coverage from test BAM file matches", {
    bam.file <- system.file("extdata", "ex1.bam", package = "PureCN", 
        mustWork = TRUE)
    interval.file <- system.file("extdata", "ex1_intervals_headered.txt", 
        package = "PureCN", mustWork = TRUE)
    coverage <- calculateBamCoverageByInterval(bam.file = bam.file, 
        interval.file = interval.file)
    expect_equal(coverage$average.coverage, c(37.49301, 43.78357, 39.10000),
        tolerance = 0.01)
    expect_equal(coverage$counts, c(568, 1158,  595), tolerance = 0.01)
}) 


test_that("Coverage output is correct", {
    x <- readCoverageFile(output.file)
    expect_equal(x$average.coverage, c(20.95205, 43.78357, 21.29271),
        tolerance = 0.01)
    expect_equal(x$counts, c(610, 1158, 636), tolerance = 0.01)
    interval.file <- system.file("extdata", "example_intervals.txt", 
        package = "PureCN")
    expect_error(correctCoverageBias(x, interval.file))
})

test_that("Reading BAM in chunks works", {
    fl <- system.file("extdata", "ex1.bam", package = "Rsamtools",
                      mustWork = TRUE)
    res0 <- scanBam(fl)[[1]] # always list-of-lists
    idx <- sort(sample(length(res0[[1]]), 300))
    idx <- idx[!is.na(res0$pos[idx])]
    x <- GRanges(seqnames = res0$rname[idx],
        IRanges(start = res0$pos[idx], end = res0$pos[idx] + 20))
    x$Gene <- "."
    x$on.target <- TRUE
    x$gc_bias <- NA
    x$mappability <- NA
    x$reptiming <- NA
    f2 <- tempfile()
    suppressWarnings(PureCN:::.writeIntervals(x, f2))
    r1 <- calculateBamCoverageByInterval(fl, f2)
    r2 <- calculateBamCoverageByInterval(fl, f2, chunks = 3)
    file.remove(f2)
    expect_equal(as.character(r1), as.character(x))
    expect_equal(as.character(r2), as.character(x))
    expect_equivalent(r1$counts, r2$counts)
})    
file.remove(output.file)
