#' Hurricane dataset for the North Atlantic region (up to 2017).
#'
#' This is a dataset of the strongest hurricanes in the North Atlantic region. The dataset is a clean up version of the dataset \code{\link[HURDAT]{AL}} which is available in the HURDAT package.
#' @format A data frame with 1831 rows and 3 variables:
#' \describe{
#'   \item{Year}{Year of the tropical cyclone occure (up to 2017)}
#'   \item{Key}{Unique key identifying the tropical cyclone. Formatted like AABBCCCC where AA is Basin, BB is YearNum and CC is Year}
#'  \item{WmaxST}{Maximum Wind Speed (in knots per hour) of the strongest hurricanes in the North Atlantic region}
#' }
#' @examples
#' data(Hurricane)
#' y=Hurricane$WmaxST
#' x=Hurricane$Year
#' plot(x,y)
"Hurricane"
