% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evolution_path_analysis.R
\name{evolution_mut_signature_cellwise_corr_samplewise}
\alias{evolution_mut_signature_cellwise_corr_samplewise}
\title{evolution_mut_signature_cellwise_corr_samplewise}
\usage{
evolution_mut_signature_cellwise_corr_samplewise(
  meta_info,
  name,
  pat_info,
  sample_info,
  clone_path_file_prefix,
  feature_file_prefix,
  output_file_prefix,
  cell_cnt_thr = 10
)
}
\arguments{
\item{meta_info}{seurat_obj@meta.data}

\item{name}{title for your data}

\item{pat_info}{variable name for patient in the "meta_info"}

\item{sample_info}{variable name for sample in the "meta_info"}

\item{clone_path_file_prefix}{prefix of clone path file}

\item{feature_file_prefix}{prefix of feature file}

\item{output_file_prefix}{prefix of outputfile}

\item{cell_cnt_thr}{threshold for excluding clones less than the given cell count}
}
\value{
Correlation analysis between the same mutational features and transcriptional features in each evolution path
}
\description{
Correlation analysis between the same mutational features and transcriptional features in each evolution path
}
\examples{
#meta_info <- seurat_obj@meta.data
meta_info <- read.table("data/meta_info_lung", sep = "\t", check.names = F)
name <- "lung"
pat_info <- "pat_collapsed"
sample_info <- "sample_name"
clone_path_file_prefix <- "data/robustclone/"
feature_file_prefix <- "data/feature/"
output_file_prefix <- "./"
evolution_mut_signature_cellwise_corr_samplewise(meta_info, name, pat_info, sample_info, clone_path_file_prefix, feature_file_prefix, output_file_prefix)

}
