% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessIntervals.R
\name{preprocessIntervals}
\alias{preprocessIntervals}
\title{Preprocess intervals}
\usage{
preprocessIntervals(
  interval.file,
  reference.file,
  output.file = NULL,
  off.target = FALSE,
  average.target.width = 400,
  min.target.width = 100,
  min.off.target.width = 20000,
  average.off.target.width = 2e+05,
  off.target.padding = -500,
  mappability = NULL,
  min.mappability = c(0.6, 0.1, 0.7),
  reptiming = NULL,
  average.reptiming.width = 1e+05,
  exclude = NULL,
  off.target.seqlevels = c("targeted", "all"),
  small.targets = c("resize", "drop")
)
}
\arguments{
\item{interval.file}{File specifying the intervals. Interval is expected in
first column in format CHR:START-END.  Instead of a file, a \code{GRanges}
object can be provided. This allows the use of BED files for example. Note
that GATK interval files are 1-based (first position of the genome is 1).
Other formats like BED files are often 0-based. The \code{import} function
will automatically convert to 1-based \code{GRanges}.}

\item{reference.file}{Reference FASTA file.}

\item{output.file}{Optionally, write GC content file.}

\item{off.target}{Include off-target regions.}

\item{average.target.width}{Split large targets to approximately this size.}

\item{min.target.width}{Make sure that target regions are of at least
this specified width. See \code{small.targets}.}

\item{min.off.target.width}{Only include off-target regions of that
size}

\item{average.off.target.width}{Split off-target regions to that}

\item{off.target.padding}{Pad off-target regions.}

\item{mappability}{Annotate intervals with mappability score. Assumed on a scale
from 0 to 1, with score being 1/(number alignments). Expected as \code{GRanges}
object with first meta column being the score. Regions outside these ranges are
ignored, assuming that \code{mappability} covers the whole accessible genome.}

\item{min.mappability}{\code{double(3)} specifying the minimum mappability score
for on-target, off-target, and chrY regions in that order. The chrY regions
are only used for sex determination in \sQuote{PureCN} and are therefore
treated differently. Requires \code{mappability}.}

\item{reptiming}{Annotate intervals with replication timing score. Expected as
\code{GRanges} object with first meta column being the score.}

\item{average.reptiming.width}{Tile \code{reptiming} into bins of specified
width.}

\item{exclude}{Any target that overlaps with this \code{GRanges} object
will be excluded.}

\item{off.target.seqlevels}{Controls how to deal with chromosomes/contigs
found in the \code{reference.file} but not in the \code{interval.file}.}

\item{small.targets}{Strategy to deal with targets smaller than
\code{min.target.width}.}
}
\value{
Returns GC content by interval as \code{GRanges} object.
}
\description{
Optimize intervals for copy number calling by tiling long intervals and by
including off-target regions. Uses \code{scanFa} from the Rsamtools package
to retrieve GC content of intervals in a reference FASTA file. If provided,
will annotate intervals with mappability and replication timing scores.
}
\examples{

reference.file <- system.file("extdata", "ex2_reference.fa",
    package = "PureCN", mustWork = TRUE)
interval.file <- system.file("extdata", "ex2_intervals.txt",
    package = "PureCN", mustWork = TRUE)
bed.file <- system.file("extdata", "ex2_intervals.bed",
    package = "PureCN", mustWork = TRUE)
preprocessIntervals(interval.file, reference.file,
    output.file = "gc_file.txt")

intervals <- import(bed.file)
preprocessIntervals(intervals, reference.file,
    output.file = "gc_file.txt")

}
\references{
Talevich et al. (2016). CNVkit: Genome-Wide Copy Number
Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol.
}
\author{
Markus Riester
}
