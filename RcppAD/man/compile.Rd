\name{compile}
\alias{compile}
\title{Compile a c++ template to DLL suitable for MakeADFun.}
\usage{
  compile(file, flags = cxxflags(file), safebounds = TRUE)
}
\arguments{
  \item{file}{c++ file.}

  \item{flags}{Character with compile flags.}

  \item{safebounds}{Turn on preprocessor flag for bound
  checking?}
}
\description{
  Compile a c++ template into a shared object file. OpenMP
  flag is set if the template is detected to be parallel.
}

