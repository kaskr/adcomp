\name{Rinterface}
\alias{Rinterface}
\title{Create minimal R-code corresponding to a cpp template.}
\usage{
  Rinterface(file)
}
\arguments{
  \item{file}{cpp template file.}
}
\description{
  Create a skeleton of required R-code once the cpp
  template is ready.
}
\examples{
file <- system.file("examples/simple.cpp", package = "TMB")
Rinterface(file)
}

