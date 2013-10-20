\name{compile}
\alias{compile}
\title{Compile a c++ template to DLL suitable for MakeADFun.}
\usage{
  compile(file, flags = "", safebounds = TRUE,
    safeunload = TRUE, openmp = isParallelTemplate(file),
    ...)
}
\arguments{
  \item{file}{c++ file.}

  \item{flags}{Character with compile flags.}

  \item{safebounds}{Turn on preprocessor flag for bound
  checking?}

  \item{safeunload}{Turn on preprocessor flag for safe DLL
  unloading?}

  \item{openmp}{Turn on openmp flag? Auto detected for
  parallel templates.}

  \item{...}{Passed as Makeconf variables.}
}
\description{
  Compile a c++ template into a shared object file. OpenMP
  flag is set if the template is detected to be parallel.
}
\details{
  RcppAD relies on R's built in functionality to create
  shared libraries independent on the platform. A template
  is compiled by \code{compile("template.cpp")}, which will
  call R's makefile with appropriate preprocessor flags.
}

