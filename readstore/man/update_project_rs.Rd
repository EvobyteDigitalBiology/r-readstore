% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsclient.R
\name{update_project_rs}
\alias{update_project_rs}
\title{update_project_rs}
\usage{
update_project_rs(
  client,
  project_id,
  name,
  description,
  metadata,
  dataset_metadata_keys
)
}
\arguments{
\item{client}{ReadStore client}

\item{name}{Project ID to return}

\item{description}{Project name to return}

\item{metadata}{List of project metadata}

\item{dataset_metadata_keys}{vector of dataset metadata keys}
}
\value{
json object (list) with project
}
\description{
Update project from the ReadStore API
}
