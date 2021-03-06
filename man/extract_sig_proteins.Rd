% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_sig_proteins.R
\name{extract_sig_proteins}
\alias{extract_sig_proteins}
\title{extract_sig_proteins}
\usage{
extract_sig_proteins(
  data,
  column_split_by,
  p_threshold = 0.05,
  fc_threshold = 1,
  pi_score_theshold = NULL
)
}
\arguments{
\item{pi_score_theshold}{}
}
\description{
This function extracts the protein IDs for all proteins that have significant differences between the groups being compared.
}
