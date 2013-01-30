\name{summarize.probesets}
\alias{summarize.probesets}
\title{summarize.probesets}
\usage{
  summarize.probesets(phylogeny.info, oligo.data, method,
    level, verbose = TRUE, rm.phylotypes = NULL,
    species.matrix = NULL, version = 1)

  summarize.probesets(phylogeny.info, oligo.data, method,
    level, verbose = TRUE, rm.phylotypes = NULL,
    species.matrix = NULL, version = 1)
}
\arguments{
  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{level}{summarization level}

  \item{verbose}{print intermediate messages}

  \item{rm.phylotypes}{Phylotypes to exclude (a list with
  fields species, L1, L2)}

  \item{species.matrix}{Optional. Provide pre-calculated
  species-level summaries to speed up computation.
  Returns:}

  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{level}{summarization level}

  \item{verbose}{print intermediate messages}

  \item{rm.phylotypes}{Phylotypes to exclude (a list with
  fields species, L1, L2)}

  \item{species.matrix}{Optional. Provide pre-calculated
  species-level summaries to speed up computation.}

  \item{version}{Use current (1) or experimental (2)
  version Returns:}
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

