% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateMappingBiasVcf.R
\name{calculateMappingBiasGatk4}
\alias{calculateMappingBiasGatk4}
\title{Calculate Mapping Bias from GATK4 GenomicsDB}
\usage{
calculateMappingBiasGatk4(
  workspace,
  reference.genome,
  min.normals = 1,
  min.normals.betafit = 7,
  min.normals.assign.betafit = 3,
  min.normals.position.specific.fit = 10,
  min.median.coverage.betafit = 5,
  num.betafit.clusters = 9,
  min.betafit.rho = 1e-04,
  max.betafit.rho = 0.2,
  AF.info.field = "AF"
)
}
\arguments{
\item{workspace}{Path to the GenomicsDB created by \code{GenomicsDBImport}}

\item{reference.genome}{Reference FASTA file.}

\item{min.normals}{Minimum number of normals with heterozygous SNP for
calculating position-specific mapping bias.}

\item{min.normals.betafit}{Minimum number of normals with heterozygous SNP
fitting a beta distribution}

\item{min.normals.assign.betafit}{Minimum number of normals with
heterozygous SNPs to assign to a beta binomal fit cluster}

\item{min.normals.position.specific.fit}{Minimum normals to use
position-specific beta-binomial fits. Otherwise only clustered fits are
used.}

\item{min.median.coverage.betafit}{Minimum median coverage of normals with
heterozygous SNP for fitting a beta distribution}

\item{num.betafit.clusters}{Maximum number of beta binomial fit clusters}

\item{min.betafit.rho}{Minimum dispersion factor rho}

\item{max.betafit.rho}{Maximum dispersion factor rho}

\item{AF.info.field}{Field in the \code{workspace} that stores the allelic
fraction}
}
\value{
A \code{GRanges} object with mapping bias and number of normal
samples with this variant.
}
\description{
Function calculate mapping bias for each variant in the provided
panel of normals GenomicsDB.
}
\examples{

\dontrun{
resources_file <- system.file("extdata", "gatk4_pon_db.tgz",
    package = "PureCN")
tmp_dir <- tempdir()
untar(resources_file, exdir = tmp_dir)
workspace <- file.path(tmp_dir, "gatk4_pon_db")
bias <- calculateMappingBiasGatk4(workspace, "hg19")
saveRDS(bias, "mapping_bias.rds")
unlink(tmp_dir, recursive=TRUE)
}

}
\author{
Markus Riester
}
