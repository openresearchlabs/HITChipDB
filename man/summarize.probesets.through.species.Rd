\name{summarize.probesets.through.species}
\alias{summarize.probesets.through.species}
\title{Description: summarize.probesets.through.species}
\usage{
  summarize.probesets.through.species(level,
    phylogeny.info, oligo.data, method, rm.phylotypes)
}
\arguments{
  \item{level}{summarization level}

  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{rm.phylotypes}{Phylotypes to exclude (a list with
  fields species, L1, L2)

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

