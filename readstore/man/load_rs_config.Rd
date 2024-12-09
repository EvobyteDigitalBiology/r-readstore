% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsconfig.R
\name{load_rs_config}
\alias{load_rs_config}
\title{load_rs_config}
\usage{
load_rs_config(
  filename = NULL,
  default_endpoint_url = NULL,
  default_fastq_extensions = NULL,
  default_output = NULL
)
}
\arguments{
\item{filename}{Path to configuration file}

\item{default_endpoint_url}{Default endpoint URL to set if not found in file or ENV}

\item{default_fastq_extensions}{Default fastq extensions to set if not found in file or ENV}

\item{default_output}{Default output format to set if not found in file or ENV}
}
\value{
List containing username, token, endpoint URL, fastq extensions and output format
}
\description{
Load ReadStore configuration from file or ENV variables
}
