\name{NEWS}
\title{News for Package \pkg{anndataR}}

\section{Version 1.2.0, Bioconductor 3.23 Release (April 2026)}{

  \subsection{New features}{
    \itemize{
      \item{
        Allow manually setting chunk size for HDF5 writes
      }
      \item{
        Add auto-chunking for HDF5 writes to improve performance
      }
      \item{
        Optimise sparse matrix reading performance by avoiding
        \code{Matrix::sparseMatrix} and constructing objects manually
      }
    }
  }

  \subsection{Bug fixes}{
    \itemize{
      \item{
        Handle unnamed \code{SingleCellExperiment} assays in \code{as_AnnData()}
        by automatically assigning names with a warning
      }
      \item{
        Improve warnings when items fail to convert between objects and add new
        warnings and checks
      }
    }
  }

  \subsection{Documentation}{
    \itemize{
      \item{
        Updates to vignettes
      }
      \item{
        Fix broken links
      }
      \item{
        Update author details
      }
    }
  }

  \subsection{Testing}{
    \itemize{
      \item{
        Add a benchmarking system to monitor function performance
      }
      \item{
        Improvements to roundtrip tests
      }
      \item{
        Improvements to HDF5 tests
      }
      \item{
        Improvements to GitHub Actions CI
      }
    }
  }
}

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
