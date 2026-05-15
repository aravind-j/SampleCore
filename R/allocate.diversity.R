#' Allocation of Entries to be Selected from Clusters/Groups based on Diversity
#' Index Estimates for Core Collection Development
#'
#' Estimate the number of entries to be allocated from each cluster/group in the
#' entire collection to construct a core collection on the basis of different
#' metrics computed from within cluster/group diversity index estimates. The
#' following strategies are implemented. \loadmathjax
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
#' @template qual-arg
#' @param method The allocation method. Either \code{"div"} for constant or
#'   \code{"div.prop"} for proportional or \code{"div.log"} for logarithmic or
#'   \code{"div.sqrt"} for square root allocation.
#' @template div-arg
#'
#' @template alloc-returns
#'
#' @importFrom DiversityStats simpson shannon mcintosh_diversity
#' @export
#'
#' @references
#'
#' \insertAllCited
#'
#' @seealso \code{\link[SampleCore]{allocate.basic}},
#'   \code{\link[SampleCore]{allocate.distance}}
#'
#' @examples
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Prepare example data
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' # Get data
#' data("cassava_EC_gp")
#' data <- cbind(genotypes = rownames(cassava_EC_gp), cassava_EC_gp)
#' row.names(data) <- NULL
#'
#' # Column names of traits
#' quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW",
#'            "AVPW", "ARSR", "SRDM")
#' qual <- c("CUAL", "LNGS", "PTLC", "DSTA", "LFRT", "LBTEF", "CBTR", "NMLB",
#'           "ANGB", "CUAL9M", "LVC9M", "TNPR9M", "PL9M", "STRP", "STRC",
#'           "PSTR")
#'
#' # Convert qualitative data columns to factor
#' data[, qual] <- lapply(data[, qual], as.factor)
#'
#' # Convert quantitative data columns to qualitative scores
#' quant_to_score5 <- function(x) {
#'
#'   brks <- unique( quantile(x,
#'                            probs = seq(0, 1, 0.2),
#'                            na.rm = TRUE))
#'   cut(x, breaks = brks,
#'       include.lowest = TRUE,
#'       labels = seq_len(length(brks) - 1))
#' }
#'
#' data[, quant] <- lapply(data[, quant], quant_to_score5)
#'
#' traits <- c(quant, qual)
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Custom diversity index functions
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' div_fun_brillouin <- function(x) {
#'   n <- tabulate(x)
#'   n <- n[n > 0]
#'   N <- sum(n)
#'   if (N <= 1) {
#'     return(0)
#'   }
#'   (lgamma(N + 1) - sum(lgamma(n + 1)))/N
#' }
#'
#' div_fun_margalef <- function(x) {
#'   tab <- tabulate(x)
#'   tab <- tab[tab > 0]
#'   S <- length(tab)
#'   N <- length(x)
#'   if (N <= 1) {
#'     return(0)
#'   }
#'   (S - 1)/log(N)
#' }
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Diversity allocation
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' ## Shannon-Weaver Diversity Index
#' dist_out_shannon1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                     group = "Cluster",
#'                     qualitative = traits,
#'                     method = "div",
#'                     div.index = "shannon", metric = "pooled",
#'                     size = 0.2)
#' dist_out_shannon1
#'
#' dist_out_shannon2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.index = "shannon", metric = "mean",
#'                      size = 0.2)
#' dist_out_shannon2
#'
#' ##  Gini-Simpson Index
#' dist_out_simpson1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.index = "simpson", metric = "pooled",
#'                      size = 0.2)
#' dist_out_simpson1
#'
#' dist_out_simpson2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.index = "simpson", metric = "mean",
#'                      size = 0.2)
#' dist_out_simpson2
#'
#' ## McIntosh Diversity Index
#' dist_out_mcintosh1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.index = "mcintosh", metric = "pooled",
#'                      size = 0.2)
#' dist_out_mcintosh1
#'
#' dist_out_mcintosh2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.index = "mcintosh", metric = "mean",
#'                      size = 0.2)
#' dist_out_mcintosh2
#'
#' ## Richness
#' dist_out_richness1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.index = richness, metric = "pooled",
#'                      size = 0.2)
#' dist_out_richness1
#'
#' dist_out_richness2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.fun = richness, metric = "mean",
#'                      size = 0.2)
#' dist_out_richness2
#'
#' ## Brillouin Diversity Index
#' dist_out_brillouin1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.fun = div_fun_brillouin, metric = "pooled",
#'                      size = 0.2)
#' dist_out_brillouin1
#'
#' dist_out_brillouin2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.fun = div_fun_brillouin, metric = "mean",
#'                      size = 0.2)
#' dist_out_brillouin2
#'
#' ## Margalef's richness Index
#' dist_out_margalef1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.fun = div_fun_margalef, metric = "pooled",
#'                      size = 0.2)
#' dist_out_margalef1
#'
#' dist_out_margalef2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div",
#'                      div.fun = div_fun_margalef, metric = "mean",
#'                      size = 0.2)
#' dist_out_margalef2
#'
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Diversity allocation & Proportional
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' ## Shannon-Weaver Diversity Index
#' dist_prop_out_shannon1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.index = "shannon", metric = "pooled",
#'                      size = 0.2)
#' dist_prop_out_shannon1
#'
#' dist_prop_out_shannon2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.index = "shannon", metric = "mean",
#'                      size = 0.2)
#' dist_prop_out_shannon2
#'
#' ##  Gini-Simpson Index
#' dist_prop_out_simpson1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.index = "simpson", metric = "pooled",
#'                      size = 0.2)
#' dist_prop_out_simpson1
#'
#' dist_prop_out_simpson2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.index = "simpson", metric = "mean",
#'                      size = 0.2)
#' dist_prop_out_simpson2
#'
#' ## McIntosh Diversity Index
#' dist_prop_out_mcintosh1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.index = "mcintosh", metric = "pooled",
#'                      size = 0.2)
#' dist_prop_out_mcintosh1
#'
#' dist_prop_out_mcintosh2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.index = "mcintosh", metric = "mean",
#'                      size = 0.2)
#' dist_prop_out_mcintosh2
#'
#' ## Richness
#' dist_out_richness1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.index = richness, metric = "pooled",
#'                      size = 0.2)
#' dist_out_richness1
#'
#' dist_out_richness2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.fun = richness, metric = "mean",
#'                      size = 0.2)
#' dist_out_richness2
#'
#' ## Brillouin Diversity Index
#' dist_out_brillouin1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.fun = div_fun_brillouin, metric = "pooled",
#'                      size = 0.2)
#' dist_out_brillouin1
#'
#' dist_out_brillouin2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.fun = div_fun_brillouin, metric = "mean",
#'                      size = 0.2)
#' dist_out_brillouin2
#'
#' ## Margalef's richness Index
#' dist_out_margalef1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.fun = div_fun_margalef, metric = "pooled",
#'                      size = 0.2)
#' dist_out_margalef1
#'
#' dist_out_margalef2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.prop",
#'                      div.fun = div_fun_margalef, metric = "mean",
#'                      size = 0.2)
#' dist_out_margalef2
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Diversity allocation & Logarithmic
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' ## Shannon-Weaver Diversity Index
#' dist_log_out_shannon1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.index = "shannon", metric = "pooled",
#'                      size = 0.2)
#' dist_log_out_shannon1
#'
#' dist_log_out_shannon2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.index = "shannon", metric = "mean",
#'                      size = 0.2)
#' dist_log_out_shannon2
#'
#' ##  Gini-Simpson Index
#' dist_log_out_simpson1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.index = "simpson", metric = "pooled",
#'                      size = 0.2)
#' dist_log_out_simpson1
#'
#' dist_log_out_simpson2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.index = "simpson", metric = "mean",
#'                      size = 0.2)
#' dist_log_out_simpson2
#'
#' ## McIntosh Diversity Index
#' dist_log_out_mcintosh1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.index = "mcintosh", metric = "pooled",
#'                      size = 0.2)
#' dist_log_out_mcintosh1
#'
#' dist_log_out_mcintosh2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.index = "mcintosh", metric = "mean",
#'                      size = 0.2)
#' dist_log_out_mcintosh2
#'
#' ## Richness
#' dist_out_richness1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.index = richness, metric = "pooled",
#'                      size = 0.2)
#' dist_out_richness1
#'
#' dist_out_richness2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.fun = richness, metric = "mean",
#'                      size = 0.2)
#' dist_out_richness2
#'
#' ## Brillouin Diversity Index
#' dist_out_brillouin1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.fun = div_fun_brillouin, metric = "pooled",
#'                      size = 0.2)
#' dist_out_brillouin1
#'
#' dist_out_brillouin2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.fun = div_fun_brillouin, metric = "mean",
#'                      size = 0.2)
#' dist_out_brillouin2
#'
#' ## Margalef's richness Index
#' dist_out_margalef1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.fun = div_fun_margalef, metric = "pooled",
#'                      size = 0.2)
#' dist_out_margalef1
#'
#' dist_out_margalef2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.log",
#'                      div.fun = div_fun_margalef, metric = "mean",
#'                      size = 0.2)
#' dist_out_margalef2
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Diversity allocation & Square root
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' ## Shannon-Weaver Diversity Index
#' dist_sqrt_out_shannon1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.index = "shannon", metric = "pooled",
#'                      size = 0.2)
#' dist_sqrt_out_shannon1
#'
#' dist_sqrt_out_shannon2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.index = "shannon", metric = "mean",
#'                      size = 0.2)
#' dist_sqrt_out_shannon2
#'
#' ##  Gini-Simpson Index
#' dist_sqrt_out_simpson1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.index = "simpson", metric = "pooled",
#'                      size = 0.2)
#' dist_sqrt_out_simpson1
#'
#' dist_sqrt_out_simpson2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.index = "simpson", metric = "mean",
#'                      size = 0.2)
#' dist_sqrt_out_simpson2
#'
#' ## McIntosh Diversity Index
#' dist_sqrt_out_mcintosh1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.index = "mcintosh", metric = "pooled",
#'                      size = 0.2)
#' dist_sqrt_out_mcintosh1
#'
#' dist_sqrt_out_mcintosh2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.index = "mcintosh", metric = "mean",
#'                      size = 0.2)
#' dist_sqrt_out_mcintosh2
#'
#' ## Richness
#' dist_out_richness1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.index = richness, metric = "pooled",
#'                      size = 0.2)
#' dist_out_richness1
#'
#' dist_out_richness2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.fun = richness, metric = "mean",
#'                      size = 0.2)
#' dist_out_richness2
#'
#' ## Brillouin Diversity Index
#' dist_out_brillouin1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.fun = div_fun_brillouin, metric = "pooled",
#'                      size = 0.2)
#' dist_out_brillouin1
#'
#' dist_out_brillouin2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.fun = div_fun_brillouin, metric = "mean",
#'                      size = 0.2)
#' dist_out_brillouin2
#'
#' ## Margalef's richness Index
#' dist_out_margalef1 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.fun = div_fun_margalef, metric = "pooled",
#'                      size = 0.2)
#' dist_out_margalef1
#'
#' dist_out_margalef2 <-
#'   allocate.diversity(data = data, names = "genotypes",
#'                      group = "Cluster",
#'                      qualitative = traits,
#'                      method = "div.sqrt",
#'                      div.fun = div_fun_margalef, metric = "mean",
#'                      size = 0.2)
#' dist_out_margalef2
#'
allocate.diversity <- function(data, names, group,
                               qualitative,
                               method = c("div", "div.prop",
                                          "div.sqrt", "div.log"),
                               div.index = c("richness", "shannon",
                                             "simpson", "mcintosh"),
                               shannon.base = exp(1),
                               div.fun = NULL,
                               log.base = exp(1),
                               metric = c("pooled", "mean"),
                               size) {

  # Check stats ----

  checks.sample.core(data = data, size = size,
                     names = names, group = group,
                     dist.mat = NULL,
                     qualitative = qualitative,
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

    div_res <- div.fun(as.factor(cassava_EC_gp$CUAL))
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

        if (is.null(div.fun)) {
          ## Richness ----
          if (div.index == "richness") {
            out <- length(unique(trt_x))
          }

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
        } else {
          ## Custom with div.fun ----
          out <- div.fun(trt_x)
        }

        return(out)

      })

      ## Pooled diversity ----
      if (metric == "pooled") {
        metric_out <- sum(trait_wise_div)
      }

      ## Mean/Average diversity ----
      if (metric == "mean") {
        metric_out <- mean(trait_wise_div)
      }

      return(metric_out)

    })

  # Freq estimation ----

  ## "div" ----
  if (method == "div") {
    freq <- group_dist_metric / sum(group_dist_metric)
  }

  ## "div.prop" ----
  if (method == "div.prop") {
    freq <- (group_dist_metric * gpsize) /
      sum(group_dist_metric * gpsize)
  }

  ## "div.log" ----
  if (method == "div.log") {
    log_gpsize <- log(gpsize, base = log.base)
    freq <- (group_dist_metric * log_gpsize) /
      sum(group_dist_metric * log_gpsize)
  }

  ## "div.sqrt" ----
  if (method == "div.sqrt") {
    freq <- (group_dist_metric * sqrt(gpsize)) /
      sum(group_dist_metric * sqrt(gpsize))
  }

  # Get output ---
  out <- setNames(round(freq * size * tcount), names(gpsize))

  return(out)

}
