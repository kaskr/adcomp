\name{compile}
\alias{compile}
\title{Compile a c++ template to DLL suitable for MakeADFun.}
\usage{
  compile(file, flags = cxxflags(file))
}
\arguments{
  \item{file}{c++ file.}

  \item{flags}{Character with compile flags.}
}
\description{
  Compile a c++ template into a shared object file. OpenMP
  flag is set if the template is detected to be parallel.
}

