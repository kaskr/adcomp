\name{gdbsource}
\alias{gdbsource}
\title{Source R-script through gdb to get backtrace.}
\usage{
  gdbsource(file)
}
\arguments{
  \item{file}{Your R script}
}
\value{
  Object of class \code{backtrace}
}
\description{
  Source R-script through gdb to get backtrace.
}
\details{
  This function is useful for debugging templates. If a
  script aborts e.g. due to an out-of-bound index operation
  it should be fast to locate the line that caused the
  problem by running \code{gdbsource(file)}.
}

