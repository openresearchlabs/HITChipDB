\name{summarize.probesets}
\alias{summarize.probesets}
\title{Description: summarize.probesets}
\usage{
  summarize.probesets(phylogeny.info, oligo.data, method,
    level, verbose = TRUE, rm.phylotypes = NULL,
    species.matrix = NULL, probe.parameters = list())
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
  species-level summaries to speed up computation.}

  \item{probe.parameters}{Optional. If probe.parameters are
  given, the summarization is based on these and model
  parameters are not estimated. A list. One element for
  each probeset with the following probe vectors:
  affinities, variances

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

