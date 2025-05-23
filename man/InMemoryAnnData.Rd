% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InMemoryAnnData.R
\name{InMemoryAnnData}
\alias{InMemoryAnnData}
\title{InMemoryAnnData}
\description{
Implementation of an in-memory AnnData object.

See \link{AnnData-usage} for details on creating and using \code{AnnData} objects.
}
\examples{
## complete example
ad <- AnnData(
  X = matrix(1:15, 3L, 5L),
  layers = list(
    A = matrix(5:1, 3L, 5L),
    B = matrix(letters[1:5], 3L, 5L)
  ),
  obs = data.frame(row.names = LETTERS[1:3], cell = 1:3),
  var = data.frame(row.names = letters[1:5], gene = 1:5)
)
ad

## minimum example
AnnData(
  obs = data.frame(row.names = letters[1:10]),
  var = data.frame(row.names = LETTERS[1:5])
)

}
\seealso{
\link{AnnData-usage} for details on creating and using \code{AnnData} objects

Other AnnData classes: 
\code{\link{AbstractAnnData}},
\code{\link{HDF5AnnData}}
}
\concept{AnnData classes}
\section{Super class}{
\code{\link[anndataR:AbstractAnnData]{anndataR::AbstractAnnData}} -> \code{InMemoryAnnData}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{X}}{See \link{AnnData-usage}}

\item{\code{layers}}{See \link{AnnData-usage}}

\item{\code{obs}}{See \link{AnnData-usage}}

\item{\code{var}}{See \link{AnnData-usage}}

\item{\code{obs_names}}{See \link{AnnData-usage}}

\item{\code{var_names}}{See \link{AnnData-usage}}

\item{\code{obsm}}{See \link{AnnData-usage}}

\item{\code{varm}}{See \link{AnnData-usage}}

\item{\code{obsp}}{See \link{AnnData-usage}}

\item{\code{varp}}{See \link{AnnData-usage}}

\item{\code{uns}}{See \link{AnnData-usage}}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-InMemoryAnnData-new}{\code{InMemoryAnnData$new()}}
\item \href{#method-InMemoryAnnData-clone}{\code{InMemoryAnnData$clone()}}
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
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="n_obs"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-n_obs'><code>anndataR::AbstractAnnData$n_obs()</code></a></span></li>
<li><span class="pkg-link" data-pkg="anndataR" data-topic="AbstractAnnData" data-id="n_vars"><a href='../../anndataR/html/AbstractAnnData.html#method-AbstractAnnData-n_vars'><code>anndataR::AbstractAnnData$n_vars()</code></a></span></li>
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
\if{html}{\out{<a id="method-InMemoryAnnData-new"></a>}}
\if{latex}{\out{\hypertarget{method-InMemoryAnnData-new}{}}}
\subsection{Method \code{new()}}{
Creates a new instance of an in-memory \code{AnnData} object.
Inherits from \link{AbstractAnnData}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InMemoryAnnData$new(
  X = NULL,
  obs = NULL,
  var = NULL,
  layers = NULL,
  obsm = NULL,
  varm = NULL,
  obsp = NULL,
  varp = NULL,
  uns = NULL,
  shape = NULL
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
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
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-InMemoryAnnData-clone"></a>}}
\if{latex}{\out{\hypertarget{method-InMemoryAnnData-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{InMemoryAnnData$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
