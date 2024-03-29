% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessory.R
\name{mutation_from_path}
\alias{mutation_from_path}
\title{mutation_from_path}
\usage{
mutation_from_path(
  path,
  sample,
  sample_clone_mut_list,
  feature_gn,
  root_exclude
)
}
\arguments{
\item{path}{evolution path}

\item{sample}{sample name}

\item{sample_clone_mut_list}{list object which contains mutational profile as vector for each clone}

\item{feature_gn}{gene list for a given feature
\itemize{
\item "all": calling all the mutation from a given clone
\item gene_vector: calling only the mutation within a given gene list
}}

\item{root_exclue}{whether user wants to exclude mutation profile of root state or not
\itemize{
\item True
\item False
}}
}
\value{
list object of mutation profile for each clone in a given evolution path
}
\description{
mutation profile for each clone in a given evolution path
}
\examples{
load("data/clone_mutation_info/LT_S21_clone_mut_list")
path_mut_gn_list <- mutation_from_path("LT_S21__29_23_27_25", "LT_S21", sample_clone_mut_list, "all", F)
path_mut_gn_list <- mutation_from_path("LT_S21__29_23_27_25", "LT_S21", sample_clone_mut_list, angiogenesis_gn_list, F)
}
