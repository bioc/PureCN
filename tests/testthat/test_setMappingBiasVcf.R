context("setMappingBiasVcf")

vcf.file <- system.file("extdata", "example_vcf.vcf.gz", package = "PureCN")
vcf <- readVcf(vcf.file, "hg19")

test_that("Mapping bias without normal panel matches", {
    vcf.bias <- round(info(setMappingBiasVcf(vcf))$MBB, digits = 3)
    expected <- rep(0.977, 2331)
    expect_equal(vcf.bias, expected)
    vcf <- readVcf(vcf.file, "hg19", param = ScanVcfParam(samples = "LIB-02240e4"))
    vcf.bias <- round(info(setMappingBiasVcf(vcf))$MBB, digits = 3)
    expected <- rep(1, 2331)
    expect_equal(vcf.bias, expected)
})

test_that("Precomputed mapping bias matches", {
    normal_panel <- system.file("extdata", "normalpanel.vcf.gz",
        package = "PureCN")
    set.seed(123)
    mb <- calculateMappingBiasVcf(normal_panel, genome = "hg19")
    ov <- findOverlaps(vcf, mb, select = "first")
    idx <- !is.na(ov)
    expect_equivalent(head(mb$pon.count[ov[idx]], 3), c(15, 5, 27))
    expect_equal(head(mb$bias[ov[idx]], 3), c(0.212590, 0.981208,
        1.020470), tolerance = 0.005)

    normal_panel_precomp <- tempfile(fileext = ".rds")
    saveRDS(mb, file = normal_panel_precomp)
    mb <- setMappingBiasVcf(vcf, mapping.bias.file = normal_panel_precomp)
    idx <- info(mb)$MBPON > 0
    expect_equal(head(info(mb)$MBPON[idx], 3), c(15, 5, 27))
    expect_equal(head(info(mb)$MBB[idx], 3), c(0.212590, 0.981208,
        1.020470), tolerance = 0.005)
    expect_equal(weighted.mean(c(0.981208, 1.020470, 0.912353),
        c(5, 27, 12)),
        info(mb)$MBB[230], tolerance = 0.001)
    file.remove(normal_panel_precomp)

    vcf.single.file <- system.file("extdata", "example_single.vcf.gz", package = "PureCN")
    expect_error(calculateMappingBiasVcf(vcf.single.file), "only a single sample")
    expect_error(setMappingBiasVcf(vcf, mapping.bias.file = normal_panel), "rds suffix")
})

test_that("GenomicsDB import works", {
    skip_if_not(requireNamespace("genomicsdb"), "genomicsdb required")
    skip_if_not(requireNamespace("jsonlite"), "jsonlite required")
    resources_file <- system.file("extdata", "gatk4_pon_db.tgz",
        package = "PureCN")
    tmp_dir <- tempdir()
    untar(resources_file, exdir = tmp_dir)
    workspace <- file.path(tmp_dir, "gatk4_pon_db")
    bias <- calculateMappingBiasGatk4(workspace, "hg19")
    expect_equal(2101, length(bias))
    unlink(tmp_dir, recursive = TRUE)
    # newer 4.2.5.0 GenomicsDB
    resources_file <- system.file("extdata", "gatk4_m2_test_pon_db.tgz",
        package = "PureCN")
    tmp_dir <- tempdir()
    untar(resources_file, exdir = tmp_dir)
    workspace <- file.path(tmp_dir, "gatk4_m2_test_pon_db")
    bias <- calculateMappingBiasGatk4(workspace, "hg38")
    expect_equivalent(bias$pon.count, c(3, 3, 2, 6, 2, 6, 5, 6, 8, 1, 1, 1, 3, 2, 1, 1, 3, 3, 3))
    unlink(tmp_dir, recursive = TRUE)

})


test_that("Issue 184_2 is fixed", {
    vcf.184.2 <- readVcf(system.file("extdata", "issue184_2.vcf.gz",
        package = "PureCN"))
    mb.184.2 <- readRDS(system.file("extdata", "issue184_2_mb.rds",
        package = "PureCN"))
    expect_equal(1, PureCN:::.findOverlapsCheckAlt(vcf.184.2, mb.184.2))
    expect_equal(2, PureCN:::.findOverlapsCheckAlt(vcf.184.2, rev(mb.184.2)))
    alt(vcf.184.2) <- DNAStringSetList("A")
    expect_equal(2, PureCN:::.findOverlapsCheckAlt(vcf.184.2, mb.184.2))
    expect_equal(1, PureCN:::.findOverlapsCheckAlt(vcf.184.2, rev(mb.184.2)))
    mb.184.2$ALT <- c("[T]", "[A]")
    expect_equal(2, PureCN:::.findOverlapsCheckAlt(vcf.184.2, mb.184.2))
    expect_equal(1, PureCN:::.findOverlapsCheckAlt(vcf.184.2, rev(mb.184.2)))
})    
test_that("Exceptions happen with wrong parameters", {
    normal_panel <- system.file("extdata", "normalpanel.vcf.gz",
        package = "PureCN")
    expect_error(
        calculateMappingBiasVcf(normal_panel, genome = "hg19", min.normals = 0),
        "min.normals (0) must be", fixed = TRUE
    )
    expect_error(
        calculateMappingBiasVcf(normal_panel, genome = "hg19", min.normals = 10),
        "min.normals (10) cannot be larger", fixed = TRUE
    )
    expect_error(
        calculateMappingBiasVcf(normal_panel, genome = "hg19", min.normals.assign.betafit = 10),
        "min.normals.assign.betafit (10) cannot be larger", fixed = TRUE
    )
    expect_error(
        calculateMappingBiasVcf(normal_panel, genome = "hg19", min.normals.betafit = 20),
        "min.normals.betafit (20) cannot be larger", fixed = TRUE
    )
    expect_error(
        calculateMappingBiasVcf(normal_panel, genome = "hg19", min.betafit.rho = 20),
        "min.betafit.rho"
    )
    expect_error(
        calculateMappingBiasVcf(normal_panel, genome = "hg19", max.betafit.rho = 20),
        "max.betafit.rho"
    )
    expect_error(
        calculateMappingBiasVcf(normal_panel, genome = "hg19", min.betafit.rho = 0.9),
        "min.betafit.rho (0.9", fixed = TRUE
    )
})   
