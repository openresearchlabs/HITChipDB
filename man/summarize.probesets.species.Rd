\name{summarize.probesets.species}
\alias{summarize.probesets.species}
\title{Description: Probeset summarization with various methods.}
\usage{
  summarize.probesets.species(phylogeny.info, oligo.data,
    method, verbose = TRUE)
}
\arguments{
  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{verbose}{print intermediate messages

  Returns:}
}
\value{
  List with two elements: summarized.matrix (summarized
  data matrix in log10 scale) and probe.parameters (only
  used with rpa, probe-level parameter estimates)
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

