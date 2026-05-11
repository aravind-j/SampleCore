#' @param div.index The diversity index to be used to estimate within
#'   cluster/group diversity.
#' @param div.fun A function to estimate diversity index from a factor vector of
#'   qualitative trait data.
#' @param metric The metric to be computed from the diversity index. Either
#'   \code{"pooled"} or \code{"mean"}.
#' @param shannon.base The logarithm base to be used for estimation of Shannon
#'   diversity index. Default is \code{exp(1)}.
