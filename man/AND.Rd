% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AND.R
\name{AND}
\alias{AND}
\alias{dAND}
\alias{pAND}
\alias{qAND}
\alias{rAND}
\title{Quantile-based asymmetric normal distribution}
\usage{
dAND(y, mu, phi, alpha)

pAND(q, mu, phi, alpha)

qAND(beta, mu, phi, alpha)

rAND(n, mu, phi, alpha)
}
\arguments{
\item{y, q}{These are each a vector of quantiles.}

\item{mu}{This is the location parameter \eqn{\mu}.}

\item{phi}{This is the scale parameter  \eqn{\phi}.}

\item{alpha}{This is the index parameter  \eqn{\alpha}.}

\item{beta}{This is a vector of probabilities.}

\item{n}{This is the number of observations, which must be a positive integer that has length 1.}
}
\value{
\code{\link{dAND}} provides the density, \code{\link{pAND}} provides the cumulative distribution function, \code{\link{qAND}} provides the quantile function, and \code{\link{rAND}} generates a random sample from the quantile-based asymmetric normal distribution.
}
\description{
Density, cumulative distribution function, quantile function and random sample generation
from the quantile-based asymmetric normal distribution (AND) introduced in Gijbels et al. (2019a).
}
\examples{
# Quantile-based asymmetric normal distribution (AND)
# Density
rnum<-rnorm(100)
dAND(y=rnum,mu=0,phi=1,alpha=.5)

# Distribution function
pAND(q=rnum,mu=0,phi=1,alpha=.5)

# Quantile function
beta<-c(0.25,0.5,0.75)
qAND(beta=beta,mu=0,phi=1,alpha=.5)

# random sample generation
rAND(n=100,mu=0,phi=1,alpha=.5)

}
\references{
{
 Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
}
}
\seealso{
\code{\link{dQBAD}},  \code{\link{pQBAD}},  \code{\link{qQBAD}},  \code{\link{rQBAD}}
}
