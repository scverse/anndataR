\name{NEWS}
\title{News for Package \pkg{anndataR}}

\section{Version 1.4.0, Bioconductor 3.24 Release (November 2026)}{

  \subsection{New features}{
    \itemize{
      \item{
        Add support for writing Zarr v3 stores, using (i) \code{zarr_format} argument of \code{write_zarr()}/\code{as_ZarrAnnData()} or (ii) the \code{anndataR.zarr_format} option. New stores are written as Zarr v3 by default.
      }
      \item{
        Add support for reading backed objects using \item{DelayedArray} matrices, including a \code{backed} argument for \code{HDF5AnnData} and \code{InMemoryAnnData}, and conversion of backed \code{AnnData} objects to \code{SingleCellExperiment}/\code{Seurat}.
      }
      \item{
        Add support for \code{nullable-string-array} elements in H5AD and Zarr.
      }
    }
  }

  \subsection{Minor changes}{
    \itemize{
      \item{
        Support the \code{anndata.AnnData} class name used by Python \code{anndata >= 0.13}
      }
      \item{
        Handle the \code{None} layer key used by Python \code{anndata >= 0.13} which stores \code{X} as a layer
      }
      \item{
        Keep \code{NA} in double and integer matrices written to H5AD instead of replacing it with \code{NaN}, so that \code{NA} and \code{NaN} can still be told apart when reading the file back in R. Python \code{anndata} reads both as \code{NaN}.
      }
      \item{
        Use zstd compression by default in \code{write_zarr()}
      }
    }
  }

  \subsection{Bug fixes}{
    \itemize{
      \item{
        Fix \code{create_zarr()} so that Zarr stores can be created at relative paths and at paths containing regex metacharacters.
      }
      \item{
        Use the format of an existing Zarr store when writing to it and truncate the store when \code{mode = "w"}, so that a store can no longer end up with a mix of Zarr v2 and v3 nodes.
      }
      \item{
        Write strings as VLen-UTF8 rather than as fixed-width, so that string arrays read back as variable length strings in Python \code{anndata}.
      }
      \item{
        Remove the workaround for LZF compression of VLen strings in HDF5.
        This is fixed in later versions of \pkg{rhdf5}.
      }
      \item{
        Fix \code{HDF5AnnData} writing logical matrices as flat 1D datasets instead of preserving their shape.
      }
      \item{
        Fix writing \code{data.frame}s in \code{obsm}/\code{varm} so that the stored index is the \code{obs_names}/\code{var_names} of the parent object instead of R's automatic row numbers
      }
      \item{
        Fix writing non-ASCII strings to H5AD files by using UTF-8 for string datasets
      }
      \item{
        Fix the error and warning messages when reading or writing an element fails. The original error is now shown instead of a \pkg{cli} parsing error.
      }
      \item{
        Fix writing integer matrices containing \code{NA} to Zarr.
      }
      \item{
        Write length 1 \code{data.frame} columns and indices as arrays instead of scalars.
      }
      \item{
        Fix writing empty arrays to Zarr v3.
      }
      \item{
        Fix \code{as_AnnData()} failing with \code{"`rownames(X)`` is not as expected"} when the assay of a \code{SingleCellExperiment} mapped to \code{X} has different dimnames than the object.
      }
  }

  \subsection{Documentation}{
    \itemize{
      \item{
        Update CITATION to include the \pkg{anndataR} publication:
        Deconinck L, Zappia L, Cannoodt R, Morgan M, Virshup I, Sang-aram C, Bredikhin D, Seurinck R, Saeys Y, scverse core (2026). "anndataR improves interoperability between R and Python in single-cell transcriptomics" Bioinformatics, 42, btag288.doi:10.1093/bioinformatics/btag288 https://doi.org/10.1093/bioinformatics/btag288.
      }
      \item{
        Minor updates and simplifications to function documentation.
      }
    }
  }

  \subsection{Testing}{
    \itemize{
      \item{
        Add Bioconductor R-universe checks
      }
      \item{
        Add Zarr functionality to benchmarking
      }
      \item{
        Check release branches against the appropriate Bioconductor release.
        Run BiocCheck on \code{RELEASE_*} branches as well as on \code{devel}.
      }
      \item{
        Add a strict mode for Python tests with \code{ANNDATAR_TEST_REQUIRE_PYTHON=true}
      }
    }
  }
}

\section{Version 1.2.0, Bioconductor 3.23 Release (April 2026)}{

  \subsection{New features}{
    \itemize{
      \item{
        Add initial Zarr support for reading and writing \code{.zarr} stores
      }
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
