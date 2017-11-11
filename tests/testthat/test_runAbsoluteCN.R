context("runAbsoluteCN")

normal.coverage.file <- system.file("extdata", "example_normal.txt", 
    package = "PureCN")
normal2.coverage.file <- system.file("extdata", "example_normal2.txt", 
    package = "PureCN")
normal.coverage.files <- c(normal.coverage.file, normal2.coverage.file)
tumor.coverage.file <- system.file("extdata", "example_tumor.txt", 
    package = "PureCN")
vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
    package = "PureCN")
seg.file <- system.file("extdata", "example_seg.txt", package = "PureCN")
cosmic.vcf.file <- system.file("extdata", "example_cosmic.vcf.gz", 
    package = "PureCN")
data(purecn.example.output)

test_that("VCF is not necessary to produce output", {
    set.seed(123)
    target.weight.file <- tempfile(fileext = ".txt")
    createTargetWeights(tumor.coverage.file, normal.coverage.files,
        target.weight.file)
    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file,
        tumor.coverage.file = tumor.coverage.file,
        candidates = purecn.example.output$candidates,
        genome = "hg19",
        args.segmentation = list(target.weight.file = target.weight.file),
        test.purity = seq(0.4, 0.7, by = 0.05), min.ploidy = 1.4, 
        max.ploidy = 2.4, max.candidate.solutions = 1)
    file.remove(target.weight.file)

    tmpFile <- tempfile(fileext = ".rds")
    saveRDS(ret, tmpFile)
    createCurationFile(tmpFile)
    readCurationFile(tmpFile)
    file.remove(tmpFile)
    chrom <- ret$results[[1]]$seg$chrom
    expect_equal(chrom[!duplicated(chrom)], 1:22)
    expect_error(callLOH(ret), "runAbsoluteCN was run without a VCF file")
    expect_error(callMutationBurden(ret), 
        "runAbsoluteCN was run without a VCF file")
})

test_that("Exceptions happen with incorrect input data", {
    expect_error(runAbsoluteCN(tumor.coverage.file = tumor.coverage.file, 
        genome = "hg19"), 
        "Need a normal coverage file if log.ratio and seg.file")
    expect_error(runAbsoluteCN(tumor.coverage.file = tumor.coverage.file, 
        min.ploidy = 0, genome = "hg19"),
        "min.ploidy or max.ploidy not within expected range")
    expect_error(runAbsoluteCN(tumor.coverage.file = tumor.coverage.file, 
        max.ploidy = 0, genome = "hg19"),
    "min.ploidy or max.ploidy not within expected range")
    expect_error(runAbsoluteCN(tumor.coverage.file = tumor.coverage.file, 
        normal.coverage.file = normal.coverage.file, max.ploidy = "a", 
        genome = "hg19"), "max.ploidy")
    expect_error(runAbsoluteCN(tumor.coverage.file = tumor.coverage.file, 
        normal.coverage.file = normal.coverage.file, min.ploidy = "a", 
        genome = "hg19"), "min.ploidy")
    expect_error(runAbsoluteCN(tumor.coverage.file = tumor.coverage.file, 
        normal.coverage.file = normal.coverage.file, test.purity = "a", 
        genome = "hg19"), "test.purity")
    expect_error(runAbsoluteCN(tumor.coverage.file, tumor.coverage.file, 
        genome = "hg19"), "Tumor and normal are identical")
    expect_error(runAbsoluteCN(tumor.coverage.file = tumor.coverage.file, 
        log.ratio = head(purecn.example.output$input$log.ratio$log.ratio), 
        genome = "hg19"),
        "Length of log.ratio different from tumor coverage")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        prior.purity = 1, genome = "hg19"),
        "prior.purity must have the same")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        min.gof = 70, genome = "hg19"), "min.gof")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        prior.purity = c(0.6, 1.1), test.purity = c(0.2, 0.6), 
        genome = "hg19"),
        "prior.purity not within expected range")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        prior.purity = c(0.6, 0.9), test.purity = c(0.2, 0.6), 
        genome = "hg19"),
        "prior.purity must add to 1. Sum is 1.5")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        genome = "hg19", max.homozygous.loss = 1.1),
        "max.homozygous.loss not within expected range")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        genome = "hg19", prior.K = -1.1),
        "prior.K not within expected range")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        genome = "hg19", prior.contamination = c(0.1, -1.1)),
        "prior.contamination not within expected range")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        genome = "hg19", iterations = 1),
        "Iterations not in the expected range from")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        genome = "hg19", iterations = 3000),
        "Iterations not in the expected range from")
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        genome = "hg19", model.homozygous = NULL), "model.homozygous")
    normalCov <- readCoverageFile(normal.coverage.file)
    expect_error(runAbsoluteCN(normalCov[sample(length(normalCov)), 
        ], tumor.coverage.file, genome = "hg19"),
        "Interval files in normal and tumor different")
})


test_that("Example data with VCF produces expected output", {
    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        tumor.coverage.file = tumor.coverage.file, vcf.file = vcf.file, 
        genome = "hg19", test.purity = seq(0.3, 0.7, by = 0.05), 
        max.candidate.solutions = 1)
    expect_equal(ret$results[[1]]$fraction.balanced, 0.21, tolerance = 0.01)
    s <- predictSomatic(ret)
    expect_equal(s$AR/s$MAPPING.BIAS, s$AR.ADJUSTED)
    expect_equal(0.65, ret$results[[1]]$purity)
    expect_equal(0.93, ret$results[[1]]$GoF, tolerance = 0.01)
    expect_equal(as.numeric(colnames(ret$candidates$all)), c(seq(0.3, 
        0.34, by = 1/50), seq(0.35, 0.7, by = 1/30)))
})

test_that("Mapping bias function works", {
    vcf <- readVcf(vcf.file, "hg19", param = ScanVcfParam(samples = "LIB-02240e4"))
    myMappingBiasTestFun <- function(vcf, ...) {
        tmp <- rep(1, nrow(vcf))
        tmp2 <- rep(0, nrow(vcf))
        idx <- as.logical(seqnames(vcf) == "chr9" & start(vcf) == 
            35811642)
        tmp[idx] <- 0.9
        tmp2[idx] <- 2
        list(bias = tmp, pon.count = tmp2)
    }
    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        tumor.coverage.file = tumor.coverage.file, args.filterVcf = list(remove.off.target.snvs = FALSE), 
        vcf.file = vcf, genome = "hg19", test.purity = seq(0.3, 
            0.7, by = 0.05), min.ploidy = 1.4, max.ploidy = 2.4, 
        max.candidate.solutions = 1, fun.setMappingBiasVcf = myMappingBiasTestFun)
    expect_equal(0.65, ret$results[[1]]$purity)
    expect_equal(as.numeric(colnames(ret$candidates$all)), c(seq(0.3, 
        0.34, by = 1/50), seq(0.35, 0.7, by = 1/30)))
    s <- predictSomatic(ret)
    expect_equal(s$AR/s$MAPPING.BIAS, s$AR.ADJUSTED)
    expect_equal(s$MAPPING.BIAS[s$chr == "chr9" & s$start == 
        35811642], 0.9)
    expect_equal(s$pon.count[s$chr == "chr9" & s$start == 35811642], 
        2)
})

test_that("Missing Gene column in gc.gene.file is handled correctly", {
    gc_data <- read.delim(gc.gene.file, as.is = TRUE)[, -3]
    output.file <- tempfile(".txt")
    write.table(gc_data, file = output.file, row.names = FALSE, sep = "\t", 
        quote = FALSE)
    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        gc.gene.file = output.file, model.homozygous = TRUE, tumor.coverage.file = tumor.coverage.file, 
        remove.off.target.snvs = TRUE, candidates = purecn.example.output$candidates, 
        min.ploidy = 1.4, max.ploidy = 2.4, vcf.file = vcf.file, 
        genome = "hg19", test.purity = seq(0.4, 0.7, by = 0.05), 
        max.candidate.solutions = 1)
    expect_true(is.na(ret$results[[1]]$gene.calls))
    expect_error(callAlterations(ret), "requires gene-level calls")
    tumor <- readCoverageFile(tumor.coverage.file)
    normal <- readCoverageFile(normal.coverage.file)
    tumor$gc_bias <- gc_data$gc_bias
    tumor$normal.average.coverage <- normal$average.coverage
    filtered <- tumor[!overlapsAny(tumor, ret$input$log.ratio)]
    expect_equal(sum(!(filtered$average.coverage < 15 | filtered$normal.average.coverage < 
        15 | filtered$targeted.base < 5 | filtered$gc_bias < 
        0.25 | filtered$gc_bias > 0.8)),0)
    file.remove(output.file)
})

test_that("Different chromosome naming styles throw exceptions", {
    vcf <- readVcf(vcf.file, "hg19")
    seqlevelsStyle(vcf) <- "ENSEMBL"
    normCov <- readCoverageFile(normal.coverage.file)
    tumorCov <- readCoverageFile(tumor.coverage.file)
    seqlevelsStyle(normCov) <- "ENSEMBL"
    seqlevelsStyle(tumorCov) <- "ENSEMBL"
    expect_error(runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        genome = "hg19", vcf.file = vcf),
        "Different chromosome names in coverage and VCF")
    expect_error(runAbsoluteCN(normal.coverage.file = normCov, 
        tumor.coverage.file = tumorCov, gc.gene.file = gc.gene.file, 
        vcf.file = vcf, genome = "hg19", test.purity = seq(0.3, 
            0.7, by = 0.05), max.candidate.solutions = 1),
        "tumor.coverage.file and gc.gene.file do not")

    output.file1 <- tempfile(fileext = ".txt")
    output.file2 <- tempfile(fileext = ".txt")
    gc1 <- read.delim(gc.gene.file, as.is = TRUE)
    gc1[, 1] <- as.character(tumorCov)
    write.table(gc1, file = output.file1, row.names = FALSE, sep = "\t", 
        quote = FALSE)
    gc2 <- gc1
    gc2$Gene[2:4] <- paste("SCNN1D,XXXL")
    write.table(gc2, file = output.file2, row.names = FALSE, sep = "\t", 
        quote = FALSE)
    set.seed(123)
    ret <- runAbsoluteCN(normal.coverage.file = normCov, tumor.coverage.file = tumorCov, 
        gc.gene.file = output.file1, plot.cnv = FALSE, min.ploidy = 1.4, 
        max.ploidy = 2.4, vcf.file = vcf, genome = "hg19", test.purity = seq(0.4, 
            0.7, by = 0.05), max.candidate.solutions = 1)
    set.seed(123)
    ret2 <- runAbsoluteCN(normal.coverage.file = normCov, tumor.coverage.file = tumorCov, 
        gc.gene.file = output.file2, plot.cnv = FALSE, min.ploidy = 1.4, 
        max.ploidy = 2.4, vcf.file = vcf, genome = "hg19", test.purity = seq(0.4, 
            0.7, by = 0.05), max.candidate.solutions = 1)
    gpnmb <- callAlterations(ret)["GPNMB", ]
    expect_equal(as.character(gpnmb$chr), "7")
    expect_true(gpnmb$start > 2.3e+07)
    expect_true(gpnmb$end < 23400000)
    expect_true(gpnmb$C >= 6)
    expect_equal(gpnmb$type, "AMPLIFICATION")
    caAret <- callAlterations(ret, all.genes = TRUE)
    caAret2 <- callAlterations(ret2, all.genes = TRUE)
    expect_equal(nrow(caAret2), nrow(caAret) + 1)
    expect_equal(caAret2[rownames(caAret), "C"], caAret$C)
    expect_equal(caAret2[rownames(caAret), "start"], caAret$start)
    expect_equal(caAret2["XXXL", "number.targets"], 3)
    expect_equal(caAret2["XXXL", "start"], 1216606)
    expect_equal(caAret2["XXXL", "end"], 1217696)
    expect_equal(caAret2["XXXL", "C"], caAret["SCNN1D", "C"])
    expect_equal(caAret2["XXXL", "seg.mean"], caAret["SCNN1D", 
        "seg.mean"])
    plotAbs(ret, 1, type = "all")
    expect_error(plotAbs(ret, 1, type = "BAF", chr = "chr1"))
    plotAbs(ret, 1, type = "BAF", chr = "1")
    file.remove(c(output.file1, output.file2))
})

test_that("Run with example seg.file works", { 
    ret <- runAbsoluteCN(tumor.coverage.file = tumor.coverage.file, 
        seg.file = seg.file, vcf.file = vcf.file, max.candidate.solutions = 1, 
        genome = "hg19", min.ploidy = 1.4, max.ploidy = 2.4, 
        test.purity = seq(0.4, 0.7, by = 0.05), sampleid = "Sample2")
    expect_equal(0.65, ret$results[[1]]$purity)
    ret <- runAbsoluteCN(seg.file = seg.file, gc.gene.file = gc.gene.file, 
        vcf.file = vcf.file, max.candidate.solutions = 1, genome = "hg19", 
        test.purity = seq(0.4, 0.7, by = 0.05), verbose = FALSE)
    expect_equal(0.65, ret$results[[1]]$purity)
    tmp <- read.delim(seg.file, as.is = TRUE)
    colnames(tmp)[1:4] <- c("Name", "Chromosome", "Start", "End")
    output.file <- tempfile(fileext = ".seg")
    write.table(tmp, file = output.file, quote = FALSE, row.names = FALSE, 
        sep = "\t")
    expect_error(runAbsoluteCN(seg.file = output.file, gc.gene.file = gc.gene.file, 
        vcf.file = vcf.file, max.candidate.solutions = 1, genome = "hg19", 
        min.ploidy = 1.4, max.ploidy = 2.4, test.purity = seq(0.4, 
            0.7, by = 0.05), verbose = FALSE))
    file.remove(output.file)
    tmp <- read.delim(seg.file, as.is = TRUE)
    tmp2 <- tmp
    tmp2[, 1] <- "Sample2"
    tmp <- rbind(tmp, tmp2)
    output.file <- tempfile(fileext = ".seg")
    write.table(tmp, file = output.file, quote = FALSE, row.names = FALSE, 
        sep = "\t")
    expect_error(runAbsoluteCN(seg.file = output.file, gc.gene.file = gc.gene.file, 
        vcf.file = vcf.file, max.candidate.solutions = 1, genome = "hg19", 
        min.ploidy = 1.4, max.ploidy = 2.4, test.purity = seq(0.4, 
            0.7, by = 0.05), verbose = FALSE),
        "seg.file contains multiple samples and sampleid missing")
    expect_error(runAbsoluteCN(seg.file = output.file, gc.gene.file = gc.gene.file, 
        sampleid = "Sample3", vcf.file = vcf.file, max.candidate.solutions = 1, 
        genome = "hg19", test.purity = seq(0.3, 0.7, by = 0.05), 
        verbose = FALSE),
        "contains multiple samples and sampleid does not match")
    ret <- runAbsoluteCN(seg.file = output.file, gc.gene.file = gc.gene.file, 
        sampleid = "Sample1", vcf.file = vcf.file, max.candidate.solutions = 1, 
        genome = "hg19", test.purity = seq(0.3, 0.7, by = 0.05), 
        verbose = FALSE, max.ploidy = 2.5)
    file.remove(output.file)

    testSeg <- function(seg, ...) return(seg)
    res <- runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        seg.file = seg.file, fun.segmentation = testSeg, min.ploidy = 1.4, 
        max.ploidy = 2.4, test.purity = seq(0.4, 0.7, by = 0.05), 
        max.candidate.solutions = 1, genome = "hg19")
    seg <- read.delim(seg.file)
    expect_equal(nrow(res$results[[1]]$seg), nrow(seg))
    expect_equal(res$results[[1]]$seg$seg.mean, seg$seg.mean)
})

test_that("Run with provided log-ratios works", {
    log.ratio <- calculateLogRatio(readCoverageFile(normal.coverage.file), 
        readCoverageFile(tumor.coverage.file))
    ret <- runAbsoluteCN(log.ratio = log.ratio, gc.gene.file = gc.gene.file, 
        vcf.file = vcf.file, max.candidate.solutions = 1, genome = "hg19",
        min.ploidy = 1.5, max.ploidy = 2.1, plot.cnv = FALSE,
        test.purity = seq(0.4, 0.7, by = 0.05))
    expect_equal(0.65, ret$results[[1]]$purity)
    ret <- runAbsoluteCN(log.ratio = log.ratio, seg.file = seg.file, 
        gc.gene.file = gc.gene.file, vcf.file = vcf.file, max.candidate.solutions = 1, 
        min.ploidy = 1.5, max.ploidy = 2.1,
        genome = "hg19", test.purity = seq(0.3, 0.7, by = 0.05))
    expect_equal(0.65, ret$results[[1]]$purity)
})

test_that("Betabin model runs with example data", {
    vcf <- readVcf(vcf.file, "hg19")
    gt1 <- round(geno(vcf)$DP[, 1] * as.numeric(geno(vcf)$FA[, 
        1]))
    ad1 <- lapply(seq_along(gt1), function(i) as.numeric(c(geno(vcf)$DP[i, 
        1] - gt1[i], gt1[i])))
    names(ad1) <- names(geno(vcf)$DP[, 1])
    geno(vcf)$AD[, 1] <- ad1
    geno(vcf)$FA <- NULL
    geno(vcf)$DP <- NULL
    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        tumor.coverage.file = tumor.coverage.file, sampleid = "LIB-02252e4", 
        vcf.file = vcf, max.candidate.solutions = 1, genome = "hg19", 
        min.ploidy = 1.5, max.ploidy = 2.1, plot.cnv = FALSE, 
        cosmic.vcf.file = cosmic.vcf.file, model = "betabin", 
        test.purity = seq(0.3, 0.7, by = 0.01))
    expect_equal(0.65, ret$results[[1]]$purity)
})

test_that("max.ploidy works", {
    ret <- runAbsoluteCN(normal.coverage.file, tumor.coverage.file, 
        min.ploidy = 2.2, max.ploidy = 4, genome = "hg19", test.purity = seq(0.3, 
            0.7, by = 0.05), plot.cnv = FALSE, max.candidate.solutions = 1)
    expect_true(ret$results[[1]]$ploidy <= 4)
    expect_true(ret$results[[1]]$ploidy >= 2)
})

test_that("normalDB objects are used correctly", {
    expect_error(runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        tumor.coverage.file = tumor.coverage.file, genome = "hg19", 
        normalDB = vcf.file),
        "normalDB not a valid normalDB object")

    normalDB <- createNormalDatabase(normal.coverage.files)
    tmp <- normalDB
    tmp$normal.coverage.files <- NULL
    expect_error(runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        tumor.coverage.file = tumor.coverage.file, genome = "hg19", 
        normalDB = tmp), "normalDB appears to be empty")

    ret <- runAbsoluteCN(normal.coverage.file = normal.coverage.file, 
        tumor.coverage.file = tumor.coverage.file, genome = "hg19", 
        vcf.file = vcf.file, sampleid = "Sample1", gc.gene.file = gc.gene.file, 
        normalDB = normalDB, args.filterTargets = list(filter.lowhigh.gc = 0), 
        plot.cnv = FALSE, max.ploidy = 3, test.purity = seq(0.4, 
            0.7, by = 0.05), max.candidate.solutions = 1)
    tumor <- readCoverageFile(tumor.coverage.file)
    normal <- readCoverageFile(normal.coverage.file)
    idx <- overlapsAny(tumor, ret$input$log.ratio)
    cutoff <- median(normalDB$exon.median.coverage) * 0.3
    plotAbs(ret, 1, type = "volcano")
    expect_equal(sum(!(normalDB$exon.median.coverage[!idx] < 
        cutoff | width(tumor[!idx]) < 5 | tumor$average.coverage[!idx] < 
        15 | normal$average.coverage[!idx] < 15 | normalDB$fraction.missing[!idx] > 
        0.05)), 0)
    expect_true(sum(!(normalDB$exon.median.coverage[idx] < cutoff | 
        width(tumor[idx]) < 5 | tumor$average.coverage[idx] < 
        15 | normal$average.coverage[idx] < 15 | normalDB$fraction.missing[idx] > 
        0.05)) > 9000)
})

test_that("Missing on.target column in coverage data is handled gracefully", {
    tumor <- readCoverageFile(tumor.coverage.file)
    tumor[1]$on.target <- FALSE
    x <- runAbsoluteCN(normal.coverage.file, tumor, genome = "hg19", 
        min.ploidy = 1.4, max.ploidy = 2.4, test.purity = seq(0.4, 
            0.7, by = 0.05), plot.cnv = FALSE,
        max.candidate.solutions = 1, gc.gene.file = gc.gene.file)
})
