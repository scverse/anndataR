% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDF5AnnData.R
\name{HDF5AnnData}
\alias{HDF5AnnData}
\title{HDF5AnnData}
\description{
Implementation of an HDF5-backed \code{AnnData} object.

See \link{AnnData-usage} for details on creating and using \code{AnnData} objects.
}
\seealso{
\link{AnnData-usage} for details on creating and using \code{AnnData} objects

Other AnnData classes: 
\code{\link{AbstractAnnData}},
\code{\link{InMemoryAnnData}}
}
\concept{AnnData classes}
\section{Super class}{
\code{\link[anndataR:AbstractAnnData]{anndataR::AbstractAnnData}} -> \code{HDF5AnnData}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{X}}{See \link{AnnData-usage}}

\item{\code{layers}}{See \link{AnnData-usage}}

\item{\code{obsm}}{See \link{AnnData-usage}}

\item{\code{varm}}{See \link{AnnData-usage}}

\item{\code{obsp}}{See \link{AnnData-usage}}

\item{\code{varp}}{See \link{AnnData-usage}}

\item{\code{obs}}{See \link{AnnData-usage}}

\item{\code{var}}{See \link{AnnData-usage}}

\item{\code{obs_names}}{See \link{AnnData-usage}}

\item{\code{var_names}}{See \link{AnnData-usage}}

\item{\code{uns}}{See \link{AnnData-usage}}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-HDF5AnnData-new}{\code{HDF5AnnData$new()}}
\item \href{#method-HDF5AnnData-close}{\code{HDF5AnnData$close()}}
\item \href{#method-HDF5AnnData-n_obs}{\code{HDF5AnnData$n_obs()}}
\item \href{#method-HDF5AnnData-n_vars}{\code{HDF5AnnData$n_vars()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="as_HDF5AnnData"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-as_HDF5AnnData'><code>anndataR::AbstractAnnData$as_HDF5AnnData()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="as_InMemoryAnnData"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-as_InMemoryAnnData'><code>anndataR::AbstractAnnData$as_InMemoryAnnData()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="as_Seurat"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-as_Seurat'><code>anndataR::AbstractAnnData$as_Seurat()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="as_SingleCellExperiment"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-as_SingleCellExperiment'><code>anndataR::AbstractAnnData$as_SingleCellExperiment()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="layers_keys"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-layers_keys'><code>anndataR::AbstractAnnData$layers_keys()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="obs_keys"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-obs_keys'><code>anndataR::AbstractAnnData$obs_keys()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="obsm_keys"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-obsm_keys'><code>anndataR::AbstractAnnData$obsm_keys()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="obsp_keys"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-obsp_keys'><code>anndataR::AbstractAnnData$obsp_keys()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="print"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-print'><code>anndataR::AbstractAnnData$print()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="shape"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-shape'><code>anndataR::AbstractAnnData$shape()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="to_HDF5AnnData"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-to_HDF5AnnData'><code>anndataR::AbstractAnnData$to_HDF5AnnData()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="to_InMemoryAnnData"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-to_InMemoryAnnData'><code>anndataR::AbstractAnnData$to_InMemoryAnnData()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="to_Seurat"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-to_Seurat'><code>anndataR::AbstractAnnData$to_Seurat()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="to_SingleCellExperiment"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-to_SingleCellExperiment'><code>anndataR::AbstractAnnData$to_SingleCellExperiment()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="uns_keys"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-uns_keys'><code>anndataR::AbstractAnnData$uns_keys()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="var_keys"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-var_keys'><code>anndataR::AbstractAnnData$var_keys()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="varm_keys"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-varm_keys'><code>anndataR::AbstractAnnData$varm_keys()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="varp_keys"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-varp_keys'><code>anndataR::AbstractAnnData$varp_keys()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="write_h5ad"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-write_h5ad'><code>anndataR::AbstractAnnData$write_h5ad()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-HDF5AnnData-new"></a>}}
\if{latex}{\out{\hypertarget{method-HDF5AnnData-new}{}}}
\subsection{Method \code{new()}}{
Close the HDF5 file when the object is garbage collected


\code{HDF5AnnData} constructor
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDF5AnnData$new(
  file,
  X = NULL,
  obs = NULL,
  var = NULL,
  layers = NULL,
  obsm = NULL,
  varm = NULL,
  obsp = NULL,
  varp = NULL,
  uns = NULL,
  shape = NULL,
  mode = c("a", "r", "r+", "w", "w-", "x"),
  compression = c("none", "gzip", "lzf")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{file}}{The file name (character) of the \code{.h5ad} file. If this file
already exits, other arguments must be \code{NULL}.}

\item{\code{X}}{See the \code{X} slot in \link{AnnData-usage}}

\item{\code{obs}}{See the \code{obs} slot in \link{AnnData-usage}}

\item{\code{var}}{See the \code{var} slot in \link{AnnData-usage}}

\item{\code{layers}}{See the \code{layers} slot in \link{AnnData-usage}}

\item{\code{obsm}}{See the \code{obsm} slot in \link{AnnData-usage}}

\item{\code{varm}}{See the \code{varm} slot in \link{AnnData-usage}}

\item{\code{obsp}}{See the \code{obsp} slot in \link{AnnData-usage}}

\item{\code{varp}}{See the \code{varp} slot in \link{AnnData-usage}}

\item{\code{uns}}{See the \code{uns} slot in \link{AnnData-usage}}

\item{\code{shape}}{Shape tuple (e.g. \code{c(n_obs, n_vars)}). Can be provided if
both \code{X} or \code{obs} and \code{var} are not provided.}

\item{\code{mode}}{The mode to open the HDF5 file. See \code{\link[=as_HDF5AnnData]{as_HDF5AnnData()}} for
details}

\item{\code{compression}}{The compression algorithm to use. See
\code{\link[=as_HDF5AnnData]{as_HDF5AnnData()}} for details}
}
\if{html}{\out{</div>}}
}
\subsection{Details}{
The constructor creates a new HDF5 \code{AnnData} interface object. This can
either be used to either connect to an existing \code{.h5ad} file or to
create a new one. If any additional slot arguments are set an existing
file will be overwritten.
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-HDF5AnnData-close"></a>}}
\if{latex}{\out{\hypertarget{method-HDF5AnnData-close}{}}}
\subsection{Method \code{close()}}{
Close the HDF5 file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDF5AnnData$close()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-HDF5AnnData-n_obs"></a>}}
\if{latex}{\out{\hypertarget{method-HDF5AnnData-n_obs}{}}}
\subsection{Method \code{n_obs()}}{
See the \code{n_obs} field in \link{AnnData-usage}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDF5AnnData$n_obs()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-HDF5AnnData-n_vars"></a>}}
\if{latex}{\out{\hypertarget{method-HDF5AnnData-n_vars}{}}}
\subsection{Method \code{n_vars()}}{
See the \code{n_vars} field in \link{AnnData-usage}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{HDF5AnnData$n_vars()}\if{html}{\out{</div>}}
}

}
}
