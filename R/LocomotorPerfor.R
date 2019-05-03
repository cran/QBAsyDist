#' Data on locomotor performance in small and large terrestrial mammals.
#'
#' A detailed description of these data is available in Iriarte-Diaz (2002). This dataset is also used in  Gijbels et al. (2019c).
#' For \eqn{n = 142} species of mammals measurements on their body length, body mass
#' (in kg) and maximum relative running speed were recorded. The maximum relative
#' running speed measurement takes into account the body length of the mammals, and
#' was obtained by dividing the maximum speed of the mammal species by its body
#' length.
#'
#' @format A data frame with 142 rows and 2 variables:
#' \describe{
#' \item{Body_Mass}{The body mass of \eqn{n = 142} species of mammals.}
#' \item{MRRS}{The maximum relative running speed measurement takes into account the body length of the mammals, and was obtained by dividing the maximum speed of the mammal species by its body length.}
#' }
#' @references{
#'  Gijbels, I., Karim, R. and Verhasselt, A. (2019c).  Semiparametric quantile regression using quantile-based asymmetric family of densities. Manuscript.
#'
#'
#' Iriarte-Diaz, Jose (2002). Differential scaling of locomotor performance in small and large terrestrial mammals, \emph{Journal of Experimental Biology}, \bold{205}(18), 2897--2908.
#'  }
#'
#'
#' @examples
#' data(LocomotorPerfor)
#' y=log(LocomotorPerfor$MRRS)
#' x=log(LocomotorPerfor$Body_Mass)
#' plot(x,y)
"LocomotorPerfor"
