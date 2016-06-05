test_runAbsoluteCN <- function() {
    gatk.normal.file <- system.file("extdata", "example_normal.txt", 
        package = "PureCN")
    gatk.tumor.file <- system.file("extdata", "example_tumor.txt", 
        package = "PureCN")
    vcf.file <- system.file("extdata", "example_vcf.vcf", package = "PureCN")
    gc.gene.file <- system.file("extdata", "example_gc.gene.file.txt", 
        package = "PureCN")
    data(purecn.example.output)

    # run without a VCF
    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gatk.tumor.file=gatk.tumor.file, 
        candidates=purecn.example.output$candidates, 
        max.candidate.solutions=2)

    checkEqualsNumeric(ret$results[[1]]$purity, 0.65, tolerance=0.1)
    
    # chromosomes in segmentation ordered numerically, not alphabetically
    chrom <- ret$results[[1]]$seg$chrom
    checkEqualsNumeric(1:22,chrom[!duplicated(chrom)])

    # test with gc.gene.file without symbols
    gc2 <- read.delim(gc.gene.file, as.is=TRUE)[,-3]
    write.table(gc2, file="tmp.gc", row.names=FALSE, sep="\t", quote=FALSE)
    ret <-runAbsoluteCN(gatk.normal.file=gatk.normal.file, 
        gc.gene.file="tmp.gc",
        gatk.tumor.file=gatk.tumor.file, 
        candidates=purecn.example.output$candidates, 
        vcf.file=vcf.file, genome="hg19", max.candidate.solutions=2)
}    
