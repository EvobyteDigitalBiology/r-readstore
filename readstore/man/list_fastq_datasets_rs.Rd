% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsclient.R
\name{list_fastq_datasets_rs}
\alias{list_fastq_datasets_rs}
\title{list_fastq_datasets_rs}
\usage{
list_fastq_datasets_rs(client, project_id = NULL, project_name = NULL)
}
\arguments{
\item{client}{ReadStore client}

\item{project_id}{Project ID to filter datasets for}

\item{project_name}{Project name to filter datasets for}
}
\value{
(json) list of fastq datasets
}
\description{
Get list of fastq datasets from the ReadStore API
}
