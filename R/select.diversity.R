#' Selection of Entries from Clusters/Groups on the basis of Optimized Diversity
#'
#' Select entries from cluster/groups in the entire collection which form a
#' subset with the highest trait diversity according to a either pooled or mean
#' diversity index estimate. \loadmathjax
#'
#' For each cluster/group, multiple candidate subsets are sampled randomly and
#' the subset with the highest trait diversity according to either pooled or
#' mean diversity index estimate is retained. This is similar to the
#' "Maximization" or M strategy of
#' \insertCite{schoen_conservation_1993;textual}{SampleCore}.
#'
#' Entries listed as \code{always.selected} are mandatorily included in the
#' selection. Warnings are issued if requested allocation is smaller than the
#' number of always-selected entries in a cluster/group and/or when the
#' cluster/group does not contain enough remaining entries to fulfill the
#' allocation.
#'
#' @template general-arg
#' @template sel-arg
#' @template qualquant-arg
#' @template div-arg
#' @param n.iter Integer specifying the number of random candidate subsets
#'   generated per group to optimze the diversity metric.
#'
#' @returns A named list where each element contains the selected entry
#'   identifiers for a cluster/group.
#'
#' @importFrom DiversityStats simpson shannon mcintosh_diversity
#' @export
#'
#' @references
#'
#' \insertAllCited
#'
#' @examples
select.diversity <- function(data, names, group, alloc,
                             quantitative, qualitative,
                             always.selected,
                             div.index = c("shannon", "simpson", "mcintosh"),
                             shannon.base = exp(1),
                             div.fun = NULL,
                             metric = c("mean", "pooled"),
                             n.iter = 1000) {

  checks.sample.core(data = data, size = NULL,
                     names = names, group = group,
                     dist.mat = NULL,
                     quantitative = quantitative,
                     qualitative = qualitative,
                     log.base = NULL,
                     alloc = alloc,
                     always.selected = always.selected,
                     mode = "sel")

  div.index <- match.arg(div.index)
  metric <- match.arg(metric)

  # Prepare data ----

  data[, group] <- droplevels(data[, group])

  gp_memb <- data[, group]

  if (!is.null(quantitative)) {
    # Convert quantitative to qualitative
  }

  # traits <- c(qualitative, quantitative)
  traits <- c(qualitative)

  # choose diversity function once
  div_fun_internal <- switch(
    div.index,
    shannon  = function(x) DiversityStats::shannon(x, base = shannon.base),
    simpson  = DiversityStats::gini_simpson,
    mcintosh = DiversityStats::mcintosh_diversity
  )

  if (!is.null(div.fun)) {
    div_fun_internal <- div.fun
  }

  # Split data once
  grouped_data <- split(data, data[[group]])

  # Main group-wise selection ----

  out <-
    lapply(names(alloc), function(g) {

      sub_df <- grouped_data[[g]]

      if (is.null(sub_df) || nrow(sub_df) == 0) {
        warning(
          sprintf("Group %s has no rows in data.", g),
          call. = FALSE
        )
        return(character(0))
      }

      group_accns <- sub_df[[names]]

      # Always selected accessions

      fixed_accns <-
        intersect(always.selected, group_accns)

      alloc_n <- alloc[[g]]

      rem_accns <-
        setdiff(group_accns, fixed_accns)

      n_rem <- alloc_n - length(fixed_accns)


      # Handle alloc smaller than fixed set
      if (n_rem < 0) {

        warning(
          sprintf(
            paste0('Group %s: "alloc" (%d) is smaller than number of ',
                   '"always.selected" values (%d). ',
                   'Taking only "always.selected" values.'),
            g, alloc_n, length(fixed_accns)
          ),
          call. = FALSE
        )

        return(fixed_accns)
      }


      # Handle insufficient remaining pool

      avail_rem <- length(rem_accns)

      if (avail_rem < n_rem) {

        warning(
          sprintf(
            paste0('Group %s has only %d additional accessions but "alloc" ',
                   'requires %d additional entries. ',
                   'Taking all available accessions.'),
            g, avail_rem, n_rem
          ),
          call. = FALSE
        )

        n_rem <- avail_rem

      }

      # If no remaining sampling needed
      if (n_rem == 0) {

        return(fixed_accns)

      }

      # Generate candidate subsets
      candidate_subsets <-
        replicate(n.iter,
                  c(fixed_accns, sample(rem_accns, n_rem, replace = FALSE)),
                  simplify = FALSE)


      # Score all candidate subsets

      candidate_scores <-
        vapply(candidate_subsets, function(accns) {

          idx <- match(accns, group_accns)

          trait_div <-
            vapply(traits, function(trt) {

              trt_x <- sub_df[idx, trt]
              div_fun_internal(trt_x)

            }, numeric(1))

          switch(
            metric,
            mean = mean(trait_div, na.rm = FALSE),
            sum = sum(trait_div, na.rm = FALSE)
          )

        }, numeric(1))

      candidate_subsets[[which.max(candidate_scores)]]

    })

  names(out) <- names(alloc)

  return(out)

}
