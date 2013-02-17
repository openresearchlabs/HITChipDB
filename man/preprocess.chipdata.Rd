\name{preprocess.chipdata}
\alias{preprocess.chipdata}
\title{Description: Profiling preprocessing script}
\usage{
  preprocess.chipdata(dbuser, dbpwd, dbname,
    verbose = TRUE, host = NULL, port = NULL,
    use.precalculated.phylogeny = NULL,
    summarization.methods = c("nmf", "rpa", "frpa", "sum.through.species", "ave.through.species", "rpa.direct", "sum", "ave"),
    which.projects = NULL)
}
\arguments{
  \item{dbuser}{MySQL username}

  \item{dbpwd}{MySQL password}

  \item{dbname}{MySQL database name}

  \item{verbose}{monitor processing through intermediate
  messages}

  \item{host}{host; needed with FTP connections}

  \item{port}{port; needed with FTP connections}

  \item{use.precalculated.phylogeny}{use precalculated
  phylogeny?}

  \item{summarization.methods}{List summarization methods
  to be included in output}

  \item{which.projects}{Optionally specify the projects to
  extract. All samples from these projects will be
  included.

  Returns:}
}
\value{
  Preprocessed data and parameters
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

