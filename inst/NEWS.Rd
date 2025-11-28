\name{NEWS}
\title{News for Package \pkg{anndataR}}

\section{Version 1.0.0, Bioconductor 3.22 Release (October 2025)}{
  \subsection{Introducing anndataR!}{
    \itemize{
      \item{
        \pkg{anndataR} brings the AnnData object to R. It enables native reading
        and writing of H5AD files and conversion to/from
        \code{SingleCellExperiment} and \code{Seurat} objects.
      }
    }
  }
  \subsection{Features}{
    \itemize{
      \item{
        An \pkg{R6} class to work with \code{AnnData} objects in R (either
        in-memory or on-disk)
      }
      \item{
        Natively read/write H5AD files without needing a Python environment
      }
      \item{
        Convert to/from \code{SingleCellExperiment} objects
      }
      \item{
        Convert to/from \code{Seurat} objects
      }
    }
  }
}
