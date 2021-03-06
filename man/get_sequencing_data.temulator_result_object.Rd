% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{get_sequencing_data.temulator_result_object}
\alias{get_sequencing_data.temulator_result_object}
\title{Gets sequencing data from a TEMULATOR result object}
\usage{
\method{get_sequencing_data}{temulator_result_object}(x, idx = 1, ...)
}
\arguments{
\item{x}{object of class 'temulator_result_object'.}

\item{idx}{optional index of the mutation dataset to return (see output of print(x)).}

\item{...}{additional parameters passed to assign_mutation_label.}
}
\value{
tibble object
}
\description{
Gets sequencing data from a TEMULATOR result object
}
