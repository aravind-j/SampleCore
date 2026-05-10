#' Allocation of Entries to be Selected from Clusters/Groups based on Diversity
#' Index Estimates for Core Collection Development
#'
#' Estimate the number of entries to be allocated from each cluster/group in
#' the entire collection to construct a core collection on the basis of
#' different metrics computed from within cluster/group diversity index
#' estimates. The following strategies are implemented. \loadmathjax
#' \itemize{
#' \item{Diversity}
#' \item{Diversity & Proportional}
#' \item{Diversity & Logarithmic}
#' \item{Diversity & Square root}}
#'
#' @section Details:
#'
#'   The number of entries to be chosen from each cluster is estimated either on
#'   the basis of diversity of entries within that cluster/group alone or in
#'   combination with the size of the cluster/group (See
#'   \strong{\code{Methods}}).
#'
#'   There are several methods proposed on the basis of diversity indices such
#'   as genetic multiplicity (G) dependent method based on the range of genetic
#'   diversity \insertCite{yonezawa_sampling_1995}{SampleCore}, H strategy based
#'   on Nei's gene diversity \insertCite{nei_analysis_1973}{SampleCore} and a
#'   method based on the pooled Shannon diversity index
#'   \insertCite{bisht_assessment_1999,mahajan_sampling_1999}{SampleCore}.
#'   Similarly, measures such as expected proportion of heterozygous loci per
#'   individual and effective number of alleles have also been employed as a
#'   diversity measure for determining sample size
#'   \insertCite{franco_sampling_2006}{SampleCore}.
#'
#'   The within-cluster/group diversity is estimated as either pooled or mean
#'   value of cluster/group-wise diversity indices. The following diversity
#'   indices are implemented in this function.
#'
#'   \itemize{
#'    \item{Shannon or Shannon-Weaver or Shannon-Wiener Diversity Index or
#'    Shannon entropy
#'    (\mjseqn{H}) \insertCite{shannon_mathematical_1949,peet_measurement_1974}{SampleCore}}
#'    \item{Simpson's Index of Diversity or Gini's Diversity Index or
#'    Gini-Simpson Index or Nei's Diversity Index or Nei's Variation Index
#'    (\mjseqn{D}) or Hurlbert’s probability of interspecific encounter
#'    (\mjseqn{PIE}) \insertCite{gini_variabilita_1912,gini_variabilita_1912-1,greenberg_measurement_1956,berger_diversity_1970,hurlbert_nonconcept_1971,nei_analysis_1973,peet_measurement_1974}{SampleCore}}
#'    \item{McIntosh Diversity Index (\mjseqn{D_{Mc}})
#'    \insertCite{mcintosh_index_1967,peet_measurement_1974}{DiversityStats}}}
#'
#' @template divmethods-section
#'
#' @template general-arg
#' @template size-arg
#' @template log-arg
#' @template qualquant-arg
#' @param method The allocation method. Either \code{"div"} for constant or
#'   \code{"div.prop"} for proportional or \code{"div.log"} for logarithmic or
#'   \code{"div.sqrt"} for square root allocation.
#' @param div.index The diversity index to be used to estimate within
#'   cluster/group diversity.
#' @param div.fun A function to estimate diversity index from a factor vector of
#'   qualitative trait data.
#' @param metric The metric to be computed from the diversity index. Either
#'   \code{"pooled"} or \code{"mean"}.
#' @param shannon.base The logarithm base to be used for estimation of Shannon
#'   diversity index. Default is \code{exp(1)}.
#'
#' @returns
#'
#' @importFrom DiversityStats simpson shannon mcintosh_diversity
#' @export
#'
#' @references
#'
#' \insertAllCited
#'
#' @examples
#'
allocate.diversity <- function(data, names, group,
                               quantitative, qualitative,
                               method = c("div", "div.prop",
                                          "div.sqrt", "div.log"),
                               div.index = c("shannon", "simpson", "mcintosh"),
                               shannon.base = exp(1),
                               div.fun = NULL,
                               log.base = exp(1),
                               metric = c("pooled", "mean"),
                               size) {

  # Check stats ----

  checks.sample.core(data = data, size = size,
                     names = names, group = group,
                     dist.mat = dist.mat,
                     quantitative = NULL,
                     qualitative = NULL,
                     log.base = log.base,
                     mode = "alloc")

  method <- match.arg(method)
  metric <- match.arg(metric)

  if (!is.null(div.fun)) {
    div.index <- NULL
  }

  if (!is.null(div.index)) {

    div.index <- match.arg(div.index)

  } else {

    if (!is.function(div.fun)) {
      stop('"div.fun" must be a function.')
    }

    div_res <- div.fun(cassava_EC_gp$CUAL)
    if (!(is.numeric(div_res) && length(div_res) == 1)) {
      stop('"div.fun" should return a numeric vector of unit length.')
    }
  }

  # Prepare data ----

  data[, group] <- droplevels(data[, group])

  gp_memb <- data[, group]

  if (method != "nclust") {
    # Convert quantitative to qualitative
  }

  # traits <- c(qualitative, quantitative)
  traits <- c(qualitative)

  # Basic group stats ----

  gps <- levels(data[, group])
  gpsize <- summary(data[, group])
  gpcount <- length(gps)

  tcount <- nrow(data)


  # Compute group-wise diversity metrics ----

    group_dist_metric <- sapply(unique(gp_memb), function(g) {
      idx <- which(gp_memb == g)

      # Get trait-wise diversity indices
      trait_wise_div <- sapply(traits, function(trt) {

        trt_x <- droplevels(data[idx, trt])

        ## Shannon ----
        if (div.index == "shannon") {
          out <- DiversityStats::shannon(trt_x, base = shannon.base)
        }

        ## Simpson ----
        if (div.index == "simpson") {
          out <- DiversityStats::gini_simpson(trt_x)
        }

        ## McIntosh ----
        if (div.index == "mcintosh") {
          out <- DiversityStats::mcintosh_diversity(trt_x)
        }

        ## Custom with div.fun ----
        if (!is.null(div.fun)) {
          out <- div.fun(trt_x)
        }

        return(out)

      })

      ## Pooled diversity ----
      if (metric == "pooled") {
        sum(trait_wise_div)
      }

      ## Mean/Average diversity ----
      if (metric == "mean") {
        mean(trait_wise_div)
      }

    })

  # Freq estimation ----

  ## "dist" ----
  if (method == "dist") {
    freq <- group_dist_metric / sum(group_dist_metric)
  }

  ## "dist.prop" ----
  if (method == "dist.prop") {
    freq <- (group_dist_metric * gpsize) /
      sum(group_dist_metric * gpsize)
  }

  ## "dist.log" ----
  if (method == "dist.log") {
    log_gpsize <- log(gpsize, base = log.base)
    freq <- (group_dist_metric * log_gpsize) /
      sum(group_dist_metric * log_gpsize)
  }

  ## "dist.sqrt" ----
  if (method == "dist.sqrt") {
    freq <- (group_dist_metric * sqrt(gpsize)) /
      sum(group_dist_metric * sqrt(gpsize))
  }

  # Get output ---
  out <- setNames(round(freq * size * tcount), names(gpsize))

  return(out)

}
