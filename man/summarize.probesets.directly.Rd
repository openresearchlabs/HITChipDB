\name{summarize.probesets.directly}
\alias{summarize.probesets.directly}
\title{Description: summarize.probesets.directly}
\usage{
  summarize.probesets.directly(level, phylogeny.info,
    oligo.data, method, rm.phylotypes,
    probe.parameters = list())
}
\arguments{
  \item{level}{summarization level}

  \item{phylogeny.info}{oligo - phylotype matching
  data.frame}

  \item{oligo.data}{preprocessed probes x samples data
  matrix in log10 domain}

  \item{method}{summarization method}

  \item{rm.phylotypes}{Phylotypes to exclude (a list with
  fields species, L1, L2)}

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

