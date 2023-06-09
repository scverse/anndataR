% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDF5-write.R
\name{write_h5ad_data_frame}
\alias{write_h5ad_data_frame}
\title{Write H5AD data frame}
\usage{
write_h5ad_data_frame(value, file, name, index = NULL, version = "0.2.0")
}
\arguments{
\item{value}{Value to write}

\item{file}{Path to a H5AD file or an open H5AD handle}

\item{name}{Name of the element within the H5AD file}

\item{index}{The index to write. Can either be a vector of length equal to
the number of rows in \code{values} or a single character string giving the name
of a column in \code{values}. If \code{NULL} then \code{rownames(value)} is used.}

\item{version}{Encoding version of the element to write}
}
\description{
Write a data frame to an H5AD file
}
