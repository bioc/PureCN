\name{callLOH}
\alias{callLOH}
\title{Get regions of LOH}
\description{This function provides detailed LOH information by region.}
\usage{callLOH(res, id = 1, arm.cutoff = 0.9)}
\arguments{
  \item{res}{Return object of the runAbsoluteCN() function.}
  \item{id}{Candidate solution to extract LOH from. id=1 will  
use the maximum likelihood solution.}
  \item{arm.cutoff}{Min fraction LOH on a chromosome arm to call 
whole arm events.}
}

\value{Returns data.frame with LOH regions.    }

\author{Markus Riester}




\examples{
data(purecn.example.output)
head(callLOH(purecn.example.output))
}