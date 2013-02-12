\name{run.profiling.script}
\alias{run.profiling.script}
\title{run.profiling.script}
\usage{
  run.profiling.script(dbuser, dbpwd, dbname,
    verbose = TRUE, host = NULL, port = NULL,
    precalculated.probe.parameters = FALSE)
}
\arguments{
  \item{dbuser}{MySQL username}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySQL database name (HITChip: "Phyloarray";
  MITChip: "Phyloarray_MIT"; PITChip old: "Phyloarray_PIT";
  PITChip new: "pitchipdb")}

  \item{verbose}{verbose}

  \item{host}{host; needed with FTP connections}

  \item{port}{port; needed with FTP connections}

  \item{precalculated.probe.parameters}{Use precalculated
  probe parameters?

  Returns:}
}
\value{
  Profiling parameters. Also writes output to the
  user-specified directory.
}
\description{
  Description: Profiling main script
}
\details{
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

