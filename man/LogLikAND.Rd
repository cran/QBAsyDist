% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AND.R
\name{LogLikAND}
\alias{LogLikAND}
\title{Log-likelihood function for the quantile-based asymmetric normal distribution.}
\usage{
LogLikAND(y, mu, phi, alpha)
}
\arguments{
\item{y}{This is a vector of quantiles.}

\item{mu}{This is the location parameter \eqn{\mu}.}

\item{phi}{This is the scale parameter  \eqn{\phi}.}

\item{alpha}{This is the index parameter  \eqn{\alpha}.}
}
\value{
\code{\link{LogLikAND}} provides the value of the Log-likelihood function of the quantile-based asymmetric normal distribution.
}
\description{
The log-likelihood function \eqn{\ell_n(\mu,\phi,\alpha)=\ln[L_n(\mu,\phi,\alpha)]} in the quantile-based asymmetric normal distribution
is presented in Gijbels et al. (2019a).
}
\examples{
# Example
y<-rnorm(100)
LogLikAND(y,mu=0,phi=1,alpha=0.5)

}
\references{
{
 Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
}
}
