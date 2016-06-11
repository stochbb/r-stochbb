\name{logLikelihood}
\alias{logLikelihood}
\title{
Evaluates the log likelihood for some given data.
}
\description{
This function evaluates the log likelihood of some given data for the specified random variable
using a piece-wise linear approximation of the PDF on the specified interval.
}
\usage{
 logLikelihood(X, min, max, N, samples)
}
\arguments{
 \item{X}{the random variable.}
 \item{min}{the lower bound.}
 \item{max}{the upper bound.}
 \item{N}{the number of bins.}
 \item{samples}{the samples.}
}
