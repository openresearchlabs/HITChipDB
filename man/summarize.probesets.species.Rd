\name{summarize.probesets.species}
\alias{summarize.probesets.species}
\title{Description: Probeset summarization with various methods.}
\usage{
  summarize.probesets.species(phylogeny.info, oligo.data,
    method, verbose = TRUE,
    rm.species = c("Victivallis vadensis"))

  summarize.probesets.species(phylogeny.info, oligo.data,
    method, verbose = TRUE,
    rm.species = c("Victivallis vadensis"))
}
\arguments{
  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{verbose}{print intermediate messages}

  \item{rm.species}{Species to exclude Returns:}

  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{verbose}{print intermediate messages}

  \item{rm.species}{Species to exclude Returns:}
}
\value{
  List with two elements: summarized.matrix (summarized
  data matrix in log10 scale) and probe.parameters (only
  used with rpa, probe-level parameter estimates)

  List with two elements: summarized.matrix (summarized
  data matrix in log10 scale) and probe.parameters (only
  used with rpa, probe-level parameter estimates)
}
\description{
  Arguments:

  Arguments:
}
\author{
  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}

  Contact: Leo Lahti
  \email{microbiome-admin@googlegroups.com}
}
\references{
  See citation("microbiome")

  See citation("microbiome")
}
\keyword{utilities}

