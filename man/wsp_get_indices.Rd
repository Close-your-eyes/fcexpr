% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wsp_get_indices.R
\name{wsp_get_indices}
\alias{wsp_get_indices}
\title{Get indices for all gated populations in a flowjo workspace}
\usage{
wsp_get_indices(
  wsp,
  FCS.file.folder = NULL,
  groups = NULL,
  invert_groups = F,
  samples = NULL,
  invert_samples = F,
  lapply_fun = lapply,
  ...
)
}
\arguments{
\item{wsp}{vector of paths to flowjo workspaces}

\item{FCS.file.folder}{path to folder(s) of FCS files; may be one path for all wsp or a vector of paths, one for each wsp;
if not provided (NULL) fcs file paths are derived individually from wsps (xml)}

\item{groups}{vector or list of groups in flowjo to consider; if a list, each index corresponds to the index in wsp;
if NULL samples from all groups are read}

\item{invert_groups}{logical whether to invert group selection}

\item{samples}{vector or list of samples to select (names of FCS files), each index corresponds to the index in wsp;
if NULL all samples (from selected groups) are read}

\item{invert_samples}{logical whether to invert sample selection}

\item{lapply_fun}{lapply function name, unquoted; lapply, pbapply::pblapply or parallel::mclapply are suggested}

\item{...}{additional argument to the lapply function; mainly mc.cores when parallel::mclapply is chosen}
}
\value{
list of matrices, one for each selected sample; each column of the matrices holds the indices for the respective population
}
\description{
Every row in a fcs file represents an event. Every column is a parameter (channel). Gates select a subset (rows) of events by applying limits to usually one or two channels.
A sub-population of gated events hence may be defined by a vector of (row-) indices (TRUE - events is included in gate; FALSE - event is not included). A whole gating tree
may be represented by a matrix with n columns for n gates and m rows for m events. The output of this function may be saved to disk and applied to fcs files with fcexpr::inds_get_ff
in order to obtain subsetted flowfframes representing gated populations in flowjo.
}
\details{
Geometric gate definitions from flowjo are applied with CytoML::flowjo_to_gatingset and indices matrices are obtained with flowWorkspace::gh_pop_get_indices_mat.
This process may take a while depend upon size of fcs files as a .h5 file is written to disk for every fcs file before indices can be derived. Hence, it is recommended
to save the indices-matrices to disk in case of large FCS files.
}
\examples{
\dontrun{
ind_mat <- fcexpr::wsp_get_indices("mypath/my.wsp")
}
}
