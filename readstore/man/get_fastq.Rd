% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readstore.R
\name{get_fastq}
\alias{get_fastq}
\title{get_fastq}
\usage{
get_fastq(client, dataset_id = NULL, dataset_name = NULL, return_type = NULL)
}
\arguments{
\item{client}{ReadStore client}

\item{dataset_id}{Dataset ID to return}

\item{dataset_name}{Dataset name to return}

\item{return_type}{The return type (currently only json)}
}
\value{
json object (list) with fastq dataset
}
\description{
Get fastq files from the ReadStore API for a dataset
ID or Name of dataset must be provided
}
