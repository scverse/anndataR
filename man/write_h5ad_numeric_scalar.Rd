% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDF5-write.R
\name{write_h5ad_numeric_scalar}
\alias{write_h5ad_numeric_scalar}
\title{Write H5AD numeric scalar}
\usage{
write_h5ad_numeric_scalar(value, file, name, version = "0.2.0")
}
\arguments{
\item{value}{Value to write}

\item{file}{Path to a H5AD file or an open H5AD handle}

\item{name}{Name of the element within the H5AD file}

\item{version}{Encoding version of the element to write}
}
\description{
Write a numeric scalar to an H5AD file
}
