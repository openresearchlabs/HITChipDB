\name{run.profiling.script}
\alias{run.profiling.script}
\title{run.profiling.script}
\usage{
  run.profiling.script(dbuser, dbpwd, dbname,
    verbose = TRUE, host = NULL, port = NULL,
    summarization.methods = c("sum", "ave", "nmf", "rpa", "frpa"),
    which.projects = NULL)
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

  \item{summarization.methods}{List summarization methods
  to be included in output}

  \item{which.projects}{Optionally specify the projects to
  extract. All samples from these projects will be
  included.

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

