% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/muscat_de.R
\name{get_FDR_empirical_plots}
\alias{get_FDR_empirical_plots}
\title{Get diagnostic plots of the empirical null.}
\usage{
get_FDR_empirical_plots(de_output_tidy, cluster_id_oi, contrast_oi)
}
\arguments{
\item{de_output_tidy}{Data frame of DE results, containing at least the following columns: cluster_id, contrast, p_val, logFC.}

\item{cluster_id_oi}{Indicate which celltype DE results should be filtered for.}

\item{contrast_oi}{Indicate which contrast DE results should be filtered for.}
}
\value{
plot object
}
\description{
\code{get_FDR_empirical_plots}  Get diagnostic plots of the empirical null.. This is the function that works under  the hood of `get_FDR_empirical_plots_all`. Credits to Jeroen Gillis (cf satuRn package)
}
