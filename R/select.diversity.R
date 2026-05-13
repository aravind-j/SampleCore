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
#'   generated per group to optimze the diversity for random search
#'   (\code{search = "random"}).
#' @param max.iter The maximum number of 1-opt passes for greedy search
#'   (\code{search = "greedy"}).
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
                             search = c("random", "greedy"),
                             n.iter = 1000,
                             max.iter = 30) {

  div.index <- match.arg(div.index)
  metric <- match.arg(metric)
  search <- match.arg(search)

  checks.sample.core(data = data, size = NULL,
                     names = names, group = group,
                     dist.mat = NULL,
                     quantitative = quantitative,
                     qualitative = qualitative,
                     log.base = NULL,
                     alloc = alloc,
                     always.selected = always.selected,
                     mode = "sel")


  # Prepare data ----

  traits <- c(qualitative)

  data[, group] <- droplevels(data[, group])

  # gp_memb <- data[, group]

  # Split data once
  grouped_data <- split(data, data[[group]])

  # Choose diversity function once ----
  div_fun_internal <-
    switch(
      div.index,
      shannon  = function(x) DiversityStats::shannon(x, base = shannon.base),
      simpson  = DiversityStats::gini_simpson,
      mcintosh = DiversityStats::mcintosh_diversity
    )

  if (!is.null(div.fun)) {
    div_fun_internal <- div.fun
  }

  # Main group-wise selection ----

  ## Random search ----

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

      ## Random search ----
      if (search ==  "random") {

        # Candidate subsets fetched randomly
        candidate_subsets <-
          replicate( n.iter, c(fixed_accns,
                               sample(rem_accns, n_rem, replace = FALSE)),
                     simplify = FALSE)

        # Get candidate subset scores
        candidate_scores <-
          vapply(candidate_subsets,
                 function(x) {
                   compute_score(accns = x,
                                 group_accns = group_accns, sub_df = sub_df,
                                 traits = traits, div_fun = div_fun_internal,
                                 metric = metric)
                 }, numeric(1))

        best_subset <- candidate_subsets[[which.max(candidate_scores)]]

      } else if (search == "greedy") { ## Greedy search ----

        ### Greedy initialization ----

        # Ignores max.iter

        selected <- fixed_accns # start from always-selected set
        pool     <- rem_accns # remaining candidates

        for (i in seq_len(n_rem)) {
          # Score each candidate added to the current selected set
          scores <- vapply(pool, function(cand) {
            compute_score(accns = c(selected, cand),
                          group_accns = group_accns,
                          sub_df = sub_df, traits = traits,
                          div_fun = div_fun_internal,
                          metric = metric)
          }, numeric(1))

          best_cand <- pool[which.max(scores)]
          selected  <- c(selected, best_cand)
          pool      <- setdiff(pool, best_cand)
        }


        ### 1-opt local search ----

        current_sel <- selected
        current_score <- compute_score(accns = selected,
                                       group_accns = group_accns,
                                       sub_df = sub_df, traits = traits,
                                       div_fun = div_fun_internal,
                                       metric = metric)

        # Accessions that are free to be swapped out
        # (i.e. not always-selected)
        swappable_sel <- setdiff(current_sel, fixed_accns)

        iter_1opt <- 0L
        repeat {

          if (iter_1opt >= max.iter) break # cap check

          iter_1opt <- iter_1opt + 1L

          best_delta <- 0  #improvement over current score
          best_swap  <- NULL # list(out_s = x, in_s = y)

          # Pool of candidates that could replace out_acc
          candidate_pool <- setdiff(rem_accns, current_sel)

          for (out_acc in swappable_sel) {
            for (in_acc in candidate_pool) {
              trial_sel <- c(setdiff(current_sel, out_acc), in_acc)
              trial_score <-
                compute_score(accns = trial_sel,
                              group_accns = group_accns,
                              sub_df = sub_df,
                              traits = traits,
                              div_fun = div_fun_internal,
                              metric = metric)
              delta <- trial_score - current_score
              if (delta > best_delta) {
                best_delta <- delta
                best_swap  <- list(out_s = out_acc, in_s = in_acc)
              }
            }
          }

          # No improving swap found – local optimum reached
          if (is.null(best_swap)) break # local optimum — natural exit

          # Apply the best swap found in this pass
          current_sel <- c(setdiff(current_sel, best_swap$out_s),
                           best_swap$in_s)
          current_score <- current_score + best_delta

          # Update swappable set after the swap
          swappable_sel <- setdiff(current_sel, fixed_accns)

        }

        best_subset <- current_sel

      }
    })

  return(out)

}


compute_score <- function(accns, group_accns, sub_df,
                          traits, div_fun, metric) {
  idx <- match(accns, group_accns)
  trait_div <- vapply(traits, function(trt) {
    div_fun(sub_df[idx, trt])
  }, numeric(1))
  switch(metric,
         mean   = mean(trait_div, na.rm = TRUE),
         pooled = sum(trait_div,  na.rm = TRUE)
  )
}




