% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ALoD.R
\name{momentALoD}
\alias{momentALoD}
\alias{meanALoD}
\alias{varALoD}
\alias{skewALoD}
\alias{kurtALoD}
\title{Moments estimation for the quantile-based asymmetric logistic distribution.}
\usage{
meanALoD(mu, phi, alpha)

varALoD(mu, phi, alpha)

skewALoD(alpha)

kurtALoD(alpha)

momentALoD(phi, alpha, r)
}
\arguments{
\item{mu}{This is the location parameter \eqn{\mu}.}

\item{phi}{This is the scale parameter  \eqn{\phi}.}

\item{alpha}{This is the index parameter  \eqn{\alpha}.}

\item{r}{This is a value which is used to calculate the \eqn{r}th moment about \eqn{\mu}.}
}
\value{
\code{\link{meanALoD}} provides the mean, \code{\link{varALoD}} provides the variance, \code{\link{skewALoD}} provides the skewness, \code{\link{kurtALoD}} provides the kurtosis, and  \code{\link{momentALoD}} provides the \eqn{r}th moment about the location parameter \eqn{\mu} of the quantile-based asymmetric logistic distribution.
}
\description{
Mean, variance, skewness, kurtosis and moments about the location parameter (i.e., \eqn{\alpha}th quantile) of the quantile-based asymmetric logistic distribution defined in Gijbels et al. (2019a) useful for quantile regression with location parameter equal to \eqn{\mu}, scale parameter \eqn{\phi} and index parameter \eqn{\alpha}.
}
\examples{
# Example
meanALoD(mu=0,phi=1,alpha=0.5)
varALoD(mu=0,phi=1,alpha=0.5)
skewALoD(alpha=0.5)
kurtALoD(alpha=0.5)
momentALoD(phi=1,alpha=0.5,r=1)


}
\references{
{
 Gijbels, I., Karim, R. and Verhasselt, A. (2019a). On quantile-based asymmetric family of distributions: properties and inference. \emph{International Statistical Review}, \url{https://doi.org/10.1111/insr.12324}.
}
}
