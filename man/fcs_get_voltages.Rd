% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_get_voltages.R
\name{fcs_get_voltages}
\alias{fcs_get_voltages}
\title{Obtain PMT-voltages from a FCS file}
\usage{
fcs_get_voltages(file_path)
}
\arguments{
\item{file_path}{path to the fcs file}
}
\value{
a data.frame with different columns depending on the machine the fcs was acquired with;
usually the 'V'-column indicates PMT voltages
}
\description{
If there is a large number of fcs files and one wants to make sure that all files
have been acquired with the same settings this functions may help by retrieving the
voltage of PMTs.
}
\examples{
\dontrun{
# When the script is saved to R_scripts in the experiment folder,
get the absolute path to the folder
wd <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
# Find workspaces
ws <- list.files(wd, pattern = '\\\\.wsp$', recursive = T, full.names = T)
# pull out paths to fcs files, optionally filter for relevant groups or so
paths <- wsx_get_fcs_paths(ws[1])
# iterate (loop) through the paths to get the voltage of each fcs file
volts <- do.call(rbind, lapply(paths, function(x) {
v <- voltage_from_fcs(x)
v$file <- basename(x)
return(v)
}))
}
}
