\name{prune16S}
\alias{prune16S}
\title{Description: filter 16S data}
\usage{
  prune16S(full16S, pmTm.margin = 2.5, complement = 1,
    mismatch = 0, rmoligos = NULL,
    remove.nonspecific.oligos = FALSE)
}
\arguments{
  \item{full16S}{full16S}

  \item{pmTm.margin}{pmTm margin}

  \item{complement}{logical}

  \item{mismatch}{logical}

  \item{rmoligos}{oligos to exclude}

  \item{remove.nonspecific.oligos}{Logical. Remove oligos
  with multiple targets.

  Returns:}
}
\value{
  filtered 16S data
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")
}
\keyword{utilities}

