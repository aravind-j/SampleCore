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
                             always.selected = NULL,
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

  SampleCore.debug <- getOption("SampleCore.debug", default = FALSE)

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

      # sub_df -> traits_mat
      traits_mat <- do.call(cbind, lapply(traits, function(tr) {
        x <- sub_df[[tr]]
        if (is.factor(x)) as.double(as.integer(x)) else as.double(x)
      }))
      colnames(traits_mat) <- traits

      ## Random search ----
      if (search ==  "random") {

        # Candidate subsets fetched randomly
        candidate_subsets <-
          replicate(n.iter, c(fixed_accns,
                               sample(rem_accns, n_rem, replace = FALSE)),
                     simplify = FALSE)

       # Get candidate subset scores
        candidate_scores <-
          vapply(candidate_subsets,
                 function(x) {
                   idx <- match(x, group_accns)
                   compute_score(idx = idx,
                                 traits_mat = traits_mat,
                                 div_fun = div_fun_internal,
                                 metric = metric)
                 }, numeric(1))

        if (SampleCore.debug) {
          print(range(candidate_scores))
        }

        best_subset <- candidate_subsets[[which.max(candidate_scores)]]

      } else if (search == "greedy") { ## Greedy search ----

        ### Greedy initialization ----

        # Ignores max.iter

        # when fixed_accns is NULL
        if (length(fixed_accns) == 0L) {
          seed_acc <- sample(rem_accns, 1L)
          selected <- seed_acc # start from always-selected set
          pool     <- setdiff(rem_accns, seed_acc) # remaining candidates
          n_to_add <- max(0L, n_rem - 1L)
        } else {
          selected <- fixed_accns
          pool     <- rem_accns
          n_to_add <- n_rem
        }

        for (i in seq_len(n_to_add)) {
          # Score each candidate added to the current selected set
          scores <- vapply(pool, function(cand) {
            idx <- match(c(selected, cand), group_accns)
            compute_score(idx = idx,
                          traits_mat = traits_mat,
                          div_fun = div_fun_internal,
                          metric = metric)
          }, numeric(1))

          best_cand <- pool[which.max(scores)]
          selected  <- c(selected, best_cand)
          pool      <- setdiff(pool, best_cand)
        }

        ### 1-opt local search ----

        current_sel <- selected
        # idx <- match(selected, group_accns)

        # integer positions into traits_mat
        idx_lookup    <- setNames(seq_len(nrow(sub_df)), group_accns)
        current_idx   <- idx_lookup[current_sel]
        fixed_idx     <- idx_lookup[fixed_accns]
        rem_idx       <- idx_lookup[rem_accns]

        current_score <-
          compute_score(idx = idx_lookup[selected],
                        traits_mat = traits_mat,
                        div_fun = div_fun_internal,
                        metric = metric)

        iter_1opt <- 0L
        repeat {

          if (iter_1opt >= max.iter) break # cap check

          iter_1opt <- iter_1opt + 1L

          # swappable and candidate pools as integer indices
          swappable_idx  <- setdiff(current_idx, fixed_idx)
          candidate_idx  <- setdiff(rem_idx, current_idx) # recomputed each pass — current_idx mutates

          if (length(swappable_idx) == 0L || length(candidate_idx) == 0L) break

          # all (out, in) pairs — integer matrix, nrow = n_pairs
          pairs <- expand.grid(out_i = swappable_idx,
                               in_i  = candidate_idx)

          # score every pair in one vapply call
          trial_scores <- vapply(seq_len(nrow(pairs)), function(k) {
            trial_idx <- c(current_idx[current_idx != pairs$out_i[k]],
                           pairs$in_i[k])
            compute_score(trial_idx, traits_mat, div_fun_internal, metric)
          }, numeric(1))

          best_k     <- which.max(trial_scores)
          best_delta <- trial_scores[best_k] - current_score

          if (is.na(best_delta) || best_delta <= 0) break   # local optimum - natural exit

          # apply best swap found in this pass
          current_idx[current_idx == pairs$out_i[best_k]] <- pairs$in_i[best_k]
          # current_score <- current_score + best_delta
          current_score <- trial_scores[best_k]

        }

        best_subset <- group_accns[current_idx]

      }
    })

  return(out)

}


compute_score <- function(idx, traits_mat,
                          div_fun, metric) {
  trait_div <- vapply(seq_len(ncol(traits_mat)), function(t) {
    div_fun(traits_mat[idx, t, drop = TRUE])
  }, numeric(1))
  switch(metric,
         mean   = mean(trait_div,  na.rm = TRUE),
         pooled = sum(trait_div,   na.rm = TRUE)
  )
}


