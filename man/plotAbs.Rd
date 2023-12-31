% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotAbs.R
\name{plotAbs}
\alias{plotAbs}
\title{Plots for analyzing PureCN solutions}
\usage{
plotAbs(
  res,
  id = 1,
  type = c("hist", "overview", "BAF", "AF", "all"),
  chr = NULL,
  germline.only = TRUE,
  show.contour = FALSE,
  purity = NULL,
  ploidy = NULL,
  alpha = TRUE,
  show.segment.means = c("SNV", "segments", "both"),
  max.mapping.bias = 0.8,
  palette.name = "Paired",
  col.snps = "#2b6391",
  col.chr.shading = "#f0f0f0",
  ...
)
}
\arguments{
\item{res}{Return object of the \code{\link{runAbsoluteCN}} function.}

\item{id}{Candidate solutions to be plotted. \code{id=1} will draw the
plot for the maximum likelihood solution.}

\item{type}{Different types of plots. \code{hist} will plot a histogram,
assigning log-ratio peaks to integer values. \code{overview} will plot all
local optima, sorted by likelihood. \code{BAF} plots
something like a B-allele frequency plot known from SNP arrays: it plots
allele frequencies of germline variants (or most likely germline when status
is not available) against copy number. \code{AF} plots observed allelic
fractions against expected (purity), maximum likelihood (optimal
multiplicity) allelic fractions. \code{all} plots types \code{BAF} and
\code{AF} for all local optima, and is useful for generating a PDF for
manual inspection.}

\item{chr}{If \code{NULL}, show all chromosomes, otherwise only the ones
specified (\code{type="BAF"} only).}

\item{germline.only}{If \code{TRUE}, show only variants most likely being
germline in BAF plot. Useful to set to \code{FALSE} (in combination with
\code{chr}) to study potential artifacts.}

\item{show.contour}{For \code{type="overview"}, display contour plot.}

\item{purity}{Display expected integer copy numbers for purity, defaults to
purity of the solution (\code{type="hist"} and \code{"AF"} only).}

\item{ploidy}{Display expected integer copy numbers for ploidy, defaults to
ploidy of the solution (\code{type="hist"} and \code{"AF"} only).}

\item{alpha}{Add transparency to the plot if VCF contains many variants
(>2000, \code{type="AF"} and \code{type="BAF"} only).}

\item{show.segment.means}{Show segment means in germline allele frequency
plot?  If \code{both}, show SNVs and segment means. If \code{SNV} show all
SNVs. Only for \code{type="AF"}.}

\item{max.mapping.bias}{Exclude variants with high mapping bias from
plotting. Note that bias is reported on an inverse scale; a variant with
mapping bias of 1 has no bias. (\code{type="AF"} and \code{type="BAF"}
only).}

\item{palette.name}{The default \code{RColorBrewer} palette.}

\item{col.snps}{The color used for germline SNPs.}

\item{col.chr.shading}{The color used for shading alternate chromosomes.}

\item{\dots}{Additonal parameters passed to the \code{plot} function.}
}
\value{
Returns \code{NULL}.
}
\description{
This function provides various plots for finding correct purity and ploidy
combinations in the results of a \code{\link{runAbsoluteCN}} call.
}
\examples{

data(purecn.example.output)
plotAbs(purecn.example.output, type="overview")
# plot details for the maximum likelihood solution (rank 1)
plotAbs(purecn.example.output, 1, type="hist")
plotAbs(purecn.example.output, 1, type="BAF")
plotAbs(purecn.example.output, 1, type = "BAF", chr="chr2")

}
\seealso{
\code{\link{runAbsoluteCN}}
}
\author{
Markus Riester
}
