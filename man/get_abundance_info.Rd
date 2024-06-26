% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_wrappers.R
\name{get_abundance_info}
\alias{get_abundance_info}
\title{get_abundance_info}
\usage{
get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells = 10, senders_oi, receivers_oi, batches = NA)
}
\arguments{
\item{sce}{SingleCellExperiment object of the scRNAseq data of interest. Contains both sender and receiver cell types.}

\item{sample_id}{Name of the meta data column that indicates from which sample/patient a cell comes from}

\item{group_id}{Name of the meta data column that indicates from which group/condition a cell comes from}

\item{celltype_id}{Name of the column in the meta data of sce that indicates the cell type of a cell.}

\item{min_cells}{Indicates the minimal number of cells that a sample should have to be considered in the DE analysis. Default: 10. See `muscat::pbDS`.}

\item{senders_oi}{Default NULL: all celltypes will be considered as senders. If you want to select specific senders_oi: you can add this here as character vector.}

\item{receivers_oi}{Default NULL: all celltypes will be considered as receivers. If you want to select specific receivers_oi: you can add this here as character vector.}

\item{batches}{NA if no batches should be corrected for. If there should be corrected for batches during DE analysis and pseudobulk expression calculation, this argument should be the name(s) of the columns in the meta data that indicate the batch(s). Should be categorical. Pseudobulk expression values will be corrected for the first element of this vector.}
}
\value{
List containing cell type abundance plots and abundance_data data frame.
}
\description{
\code{get_abundance_info} Visualize cell type abundances.
}
\examples{
\dontrun{
library(dplyr)
sample_id = "tumor"
group_id = "pEMT"
celltype_id = "celltype"
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] \%>\% unique()  
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] \%>\% unique() 
abundance_celltype_info = get_abundance_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id =  celltype_id, min_cells = 10, senders_oi = senders_oi, receivers_oi = receivers_oi)
}

}
