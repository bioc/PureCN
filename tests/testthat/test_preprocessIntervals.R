context("preprocessIntervals")

reference.file <- system.file("extdata", "ex2_reference.fa", 
    package = "PureCN", mustWork = TRUE)
interval.file <- system.file("extdata", "ex2_intervals.txt", 
    package = "PureCN", mustWork = TRUE)
bed.file <- system.file("extdata", "ex2_intervals.bed", package = "PureCN", 
    mustWork = TRUE)

output.file <- tempfile(fileext = ".txt")

test_that("GC-bias of example reference and intervals (GATK format) matches", {
    gc <- preprocessIntervals(interval.file, reference.file, 
        output.file = output.file, min.target.width = 10)
    x <- readIntervalFile(output.file)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance = 0.001)
    expect_equal(gc$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance = 0.001)
})

test_that("GC-bias of example reference and intervals (BED format) matches", {
    x <- readIntervalFile(output.file)
    intervals <- import(bed.file)
    output.file2 <- tempfile(fileext = ".txt")
    y <- preprocessIntervals(intervals, reference.file, 
        output.file = output.file2, min.target.width = 10)
    expect_equal(y$gc_bias, x$gc_bias)
    expect_equal(as.character(y), as.character(x))
    y <- readIntervalFile(output.file2)
    file.remove(output.file2)
    expect_equal(y$gc_bias, x$gc_bias)
})

test_that("exclude option works", {
    x <- readIntervalFile(output.file)
    intervals <- import(bed.file)
    y <- preprocessIntervals(intervals, reference.file,
        exclude = intervals[2], min.target.width = 10)
    expect_equal(y$gc_bias, x$gc_bias[-2])
    expect_equal(as.character(y), as.character(x)[-2])
    expect_equal(y$gc_bias, x$gc_bias[-2])
})

test_that("Exceptions happen with wrong input", {
    expect_error(preprocessIntervals(interval.file, 
        reference.file, off.target =TRUE, off.target.padding = 5), 
        "must be negative")

    interval.file2 <- tempfile(fileext = ".txt")
    idata <- read.delim(interval.file, as.is = TRUE)
    idata[3, 1] <- "seq2:0-149"
    write.table(idata, file = interval.file2, row.names = FALSE, 
        quote = FALSE)
    expect_error(preprocessIntervals(interval.file2, 
        reference.file),
        "Interval coordinates should start at 1, not at 0")
    expect_error(preprocessIntervals(interval.file, reference.file, 
        off.target = TRUE),
        "after filtering for mappability")

    file.remove(interval.file2)
})

test_that("reptiming annotated correctly", {
    reptiming.file <- system.file("extdata", "ex2_reptiming.bed",
        package = "PureCN", mustWork = TRUE)

    reptiming <- import(reptiming.file)
    intervals <- import(bed.file)
    gr <- preprocessIntervals(intervals, reference.file, 
        reptiming = reptiming, average.reptiming.width = NULL)
    expect_equal(c(17.5, 11.0, 50, 10.0, 10.0), gr$reptiming)
    gr <- preprocessIntervals(intervals, reference.file, 
        reptiming = reptiming, average.reptiming.width = 200)
    expect_equal(c((12.5+15)/2, (15+12+0)/3, 65, 10.0, 10.0), gr$reptiming)
})

test_that("long targets are split correctly", {
    gr <- preprocessIntervals(interval.file, reference.file, 
        average.target.width = 200)
    expect_equal(c(175,175), width(gr)[2:3])
    gr <- preprocessIntervals(interval.file, reference.file)
    expect_equal(sum(c(175,175)), width(gr)[2])
})
    
test_that("Offtarget settings work as expected", {
    gc <- preprocessIntervals(interval.file, reference.file, min.target.width = 10,
        off.target = TRUE, min.off.target.width = 2, off.target.padding = -2)
    expect_equal(length(gc), 11)
    intervals <- import(bed.file)
    gc2 <- preprocessIntervals(gc, reference.file, min.target.width = 10)
    expect_equal(start(gc2), start(gc))
    expect_equal(end(gc2), end(gc))
    expect_equal(gc2$mappability, gc$mappability)
    expect_equal(gc2$gc_bias, gc$gc_bias)
    if (.Platform$OS.type != "windows") {
        mappability.file <- system.file("extdata", "ex2_mappability.bigWig", 
            package = "PureCN", mustWork = TRUE)
        mappability <- import(mappability.file)
        gcMap <- preprocessIntervals(intervals, reference.file, 
            mappability = mappability, min.target.width = 10)
        expect_equal(gcMap$mappability, c(1, 1, 0.7, 1, 1), tolerance = 0.001)
    }

    mappability.file <- system.file("extdata", "ex2_mappability.bed", 
        package = "PureCN", mustWork = TRUE)
    mappability <- import(mappability.file)
    gcMap <- preprocessIntervals(intervals, reference.file, 
        mappability = mappability)
    expect_equal(gcMap$mappability, c(1, 1, 0.7, 1, 1), tolerance = 0.001)
    expect_equal(100, min(width(gcMap)), tolerance = 0.001)

    mappability4 <- mappability
    mappability4$name <- as.character(mappability$score)
    mappability4$score <- NULL
    expect_output( gcMap <- preprocessIntervals(intervals, reference.file, 
                                            mappability = mappability4), 
                  "not numeric")
    expect_equal(gcMap$mappability, c(1, 1, 0.7, 1, 1), tolerance = 0.001)

    expect_output(gcMap <- preprocessIntervals(intervals, reference.file,
        mappability = mappability, min.mappability = c(1,1,1)),
        "Removing 1 intervals with low mappability score")
    expect_equal(gcMap$mappability, c(1, 1, 1, 1), tolerance = 0.001)
    
    # check for targets that don't overlap with mappability
    mappability2 <- mappability[-1]
    gcMap <- preprocessIntervals(intervals, reference.file, 
        mappability=mappability2)
    expect_equal(gcMap$mappability, c(1, 0.7, 1, 1), tolerance = 0.001)


    gcMap <- preprocessIntervals(intervals, reference.file, 
        mappability = mappability, off.target = TRUE, 
        off.target.padding = -5, min.off.target.width = 10)
    expect_equal(gcMap$mappability, c(1, 1, 1, 1, 0.7, 1, 1, 1, 1), 
        tolerance = 0.001)

    expect_output(gr <-preprocessIntervals(intervals[1:2], reference.file, 
        off.target = TRUE, off.target.padding = -5, 
        average.off.target.width = 100, min.off.target.width = 10), 
        "contigs from off-target regions: seq2")
    expect_equal("seq1", seqlevelsInUse(gr))

    reference.file <- system.file("extdata", "ex3_reference.fa", 
        package = "PureCN", mustWork = TRUE)
    bed.file3 <- system.file("extdata", "ex3_intervals.bed", 
        package = "PureCN", mustWork = TRUE)
    intervals3 <- import(bed.file3)
    x <- preprocessIntervals(intervals3, reference.file, min.target.width = 10)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333,
        0.48, 0.36), tolerance = 0.001)
    seqlevelsStyle(intervals3) <- "Ensembl"
    x <- preprocessIntervals(intervals3, reference.file, min.target.width = 10)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance = 0.001)
    expect_error(preprocessIntervals(intervals, reference.file),
        "Chromosome naming style of interval file")
    mappability.file3 <- system.file("extdata", "ex3_mappability.bed", 
        package = "PureCN", mustWork = TRUE)
    mappability3 <- import(mappability.file3)
    seqlevelsStyle(mappability3) <- "Ensembl"
    x <- preprocessIntervals(intervals3, reference.file, 
        mappability = mappability3, min.target.width = 10)
    expect_equal(x$gc_bias, c(0.4533333, 0.5057143, 0.5733333, 
        0.48, 0.36), tolerance = 0.001)
    expect_equal(x$mappability, c(1, 1, 0.7, 1, 1), tolerance = 0.001)
})

test_that("Small targets are dealt with properly", {
    gc <- preprocessIntervals(interval.file, reference.file, 
        min.target.width = 60)
    gc <- preprocessIntervals(interval.file, reference.file, 
        min.target.width = 60, small.targets = "drop")
    expect_equal(c(101, 301,   1, 401), start(gc))
    gc <- preprocessIntervals(interval.file, reference.file, 
        min.target.width = 200)
    expect_false(any(start(gc)) < 1)
})

file.remove(output.file)
