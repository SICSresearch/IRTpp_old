% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{model.transform}
\alias{model.transform}
\title{Model transformations
Implements parameter transforms from one parameter to others dependent on the model
Currently only 3PL supported with passings between b <-> d and c <-> cp'}
\usage{
model.transform(z, model, src, target, byrow = F)
}
\arguments{
\item{z}{The parameter matrix or named list.}

\item{model}{The model to transform}

\item{src}{The original parameter to transform}

\item{target}{The target parameter of the transform}
}
\value{
z The parameters with the transform
}
\description{
Model transformations
Implements parameter transforms from one parameter to others dependent on the model
Currently only 3PL supported with passings between b <-> d and c <-> cp'
}

