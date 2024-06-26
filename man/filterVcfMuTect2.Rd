% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterVcfMuTect2.R
\name{filterVcfMuTect2}
\alias{filterVcfMuTect2}
\title{Filter VCF MuTect2}
\usage{
filterVcfMuTect2(
  vcf,
  tumor.id.in.vcf = NULL,
  ignore = c("clustered_events", "t_lod", "str_contraction", "read_position", "position",
    "fragment_length", "multiallelic", "clipping", "strand_artifact", "strand_bias",
    "slippage", "weak_evidence", "orientation", "haplotype"),
  ...
)
}
\arguments{
\item{vcf}{\code{CollapsedVCF} object, read in with the \code{readVcf}
function from the VariantAnnotation package.}

\item{tumor.id.in.vcf}{The tumor id in the VCF file, optional.}

\item{ignore}{MuTect2 flags that mark variants for exclusion.}

\item{\dots}{Additional arguments passed to \code{\link{filterVcfBasic}}.}
}
\value{
A list with elements \code{vcf}, \code{flag} and
\code{flag_comment}.  \code{vcf} contains the filtered \code{CollapsedVCF},
\code{flag} a \code{logical(1)} flag if problems were identified, further
described in \code{flag_comment}.
}
\description{
Function to remove artifacts and low confidence/quality calls from a 
GATK4/MuTect2 generated VCF file. Also applies filters defined in 
\code{filterVcfBasic}.
}
\examples{

### This function is typically only called by runAbsolute via the 
### fun.filterVcf and args.filterVcf comments.
library(VariantAnnotation)    
vcf.file <- system.file("extdata", "example.vcf.gz", package="PureCN")
vcf <- readVcf(vcf.file, "hg19")
vcf.filtered <- filterVcfMuTect(vcf)        

}
\seealso{
\code{\link{filterVcfBasic}}
}
\author{
Markus Riester
}
