#' Selection of Entries from Clusters/Groups on the basis of Optimized Diversity
#'
#' Select entries from cluster/groups in the entire collection which form a
#' subset with the highest trait diversity according to a either pooled or mean
#' diversity index estimate. \loadmathjax
#'
#' To identify subsets with highest diversity estimates, the following
#' strategies are available. These strategies are similar to the "Maximization"
#' or M strategy of \insertCite{schoen_conservation_1993;textual}{SampleCore}.
#'
#'  \subsection{Random search / Monte Carlo Method}{For each cluster/group,
#'  multiple candidate subsets are sampled randomly and the subset with the
#'  highest trait diversity according to either pooled or mean diversity index
#'  estimate is retained. The quality of the solution improves with increasing
#'  \code{n.iter} but is not guaranteed to find the global optimum
#'  \insertCite{anatoly_zhigljavsky_stochastic_2008}{SampleCore}.
#'  }
#'
#'  \subsection{Greedy search with 1-opt}{This method builds a solution
#'  incrementally by adding the accession that maximises the diversity score at
#'  each step, starting from the \code{always.selected} accessions (or a single
#'  randomly drawn accession when there are no accessions specified in
#'  \code{always.selected} present in the particular cluster/group )
#'  \insertCite{nemhauser_analysis_1978,fisher_analysis_1978,cormen_introduction_2022}{SampleCore}.
#'  The 'greedy' solution is then refined by a 1-opt local search controlled by
#'  \code{local.search} and \code{max.iter}
#'  \insertCite{lin_computer_1965}{SampleCore}. Greedy search is deterministic
#'  given a fixed \code{always.selected} set; when there are no accessions
#'  specified in \code{always.selected} present in the particular cluster/group
#'  results may vary across runs due to the random initialisation.
#'
#'  \code{local.search = "best.improvement"} scans all possible single swaps
#'  in each pass and applies the one yielding the greatest improvement before
#'  restarting. his guarantees the steepest ascent at each pass but requires
#'  evaluating all \mjseqn{k \times (n - k)} swap pairs per pass, where
#'  \mjseqn{k} is the number of swappable accessions and \mjseqn{n - k} is the
#'  size of the candidate pool
#'  \insertCite{papadimitriou_combinatorial_1998}{SampleCore}.
#'
#'  \code{local.search = "first.improvement"} applies the first swap that
#'  improves the score and immediately restarts the search. This typically
#'  requires fewer score evaluations per pass and converges faster, but may
#'  find a different local optimum than \code{"best.improvement"}
#'  \insertCite{papadimitriou_combinatorial_1998}{SampleCore}.
#'
#'  Both strategies terminate when no improving swap exists (local optimum)
#'  or when \code{max.iter} passes have been completed.
#'  }
#'
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
#' @param search Character string specifying the search strategy used to find
#'   the subset with the highest diversity score. Either \code{"random"}
#'   (default) or \code{"greedy"} (See \strong{Details}).
#' @param local.search Character string specifying the local search strategy
#'   used in the 1-opt improvement phase of the greedy search (\code{search =
#'   "greedy"}). Either \code{"best.improvement"} (default) or
#'   \code{"first.improvement"}. Ignored when \code{search = "random"}.
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
                             div.index = c("richness", "shannon",
                                           "simpson", "mcintosh"),
                             shannon.base = exp(1),
                             div.fun = NULL,
                             metric = c("mean", "pooled"),
                             search = c("random", "greedy"),
                             local.search = c("best.improvement",
                                              "first.improvement"),
                             n.iter = 1000,
                             max.iter = 30) {

  div.index <- match.arg(div.index)
  metric <- match.arg(metric)
  search <- match.arg(search)
  local.search <- match.arg(local.search)

  if (search == "random" && !missing(local.search)) {
    warning('"local.search" is ignored when search = "random"',
            call. = FALSE)
  }

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
      richness = function(x) length(unique(x)),
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

        if (SampleCore.debug) {
          message("--Greedy initialization started.---------------------\n\n")
        }

        # Ignores max.iter

        # when fixed_accns is NULL
        if (length(fixed_accns) == 0L) {
          seed_acc <- sample(rem_accns, 1L)
          selected <- seed_acc # start from always-selected set
          pool <- setdiff(rem_accns, seed_acc) # remaining candidates
          n_to_add <- max(0L, n_rem - 1L)
        } else {
          selected <- fixed_accns
          pool <- rem_accns
          n_to_add <- n_rem
        }

        idx_lookup   <- setNames(seq_len(nrow(sub_df)), group_accns)
        selected_idx <- idx_lookup[selected]
        pool_idx     <- idx_lookup[pool]

        for (i in seq_len(n_to_add)) {
          scores <- vapply(pool_idx, function(cand_i) {
            compute_score(c(selected_idx, cand_i), traits_mat, div_fun_internal, metric)
          }, numeric(1))

          best_pos     <- which.max(scores)
          selected_idx <- c(selected_idx, pool_idx[best_pos])
          pool_idx     <- pool_idx[-best_pos]   # integer remove by position — faster than setdiff
        }

        ### 1-opt local search ----

        if (SampleCore.debug) {
          message("--Local search started.------------------------------\n\n")
        }

        current_idx <- selected_idx
        fixed_idx <- idx_lookup[fixed_accns]
        rem_idx <- idx_lookup[rem_accns]
        current_score <- compute_score(current_idx, traits_mat,
                                       div_fun_internal, metric)

        # Initialize indices ONCE
        swappable_idx <- setdiff(current_idx, fixed_idx)
        candidate_idx <- setdiff(rem_idx, current_idx)

        iter_1opt <- 0L
        repeat {
          if (iter_1opt >= max.iter) break # cap check

          iter_1opt <- iter_1opt + 1L

          # Exit if no swaps are possible
          if (length(swappable_idx) == 0L || length(candidate_idx) == 0L) break

          improved  <- FALSE

          #### Best-Improvement Strategy ----
          if (local.search == "best.improvement") {
            best_overall_score <- current_score
            best_out_val <- NULL
            best_in_val  <- NULL
            best_out_pos_in_swappable <- NULL
            best_in_pos_in_candidate  <- NULL

            # Nested Loops: Scanning all possible swaps
            for (i in seq_along(swappable_idx)) {
              out_val <- swappable_idx[i]
              # Pre-calculate the subset excluding the 'out' candidate
              subset_minus_out <- current_idx[current_idx != out_val]

              for (j in seq_along(candidate_idx)) {
                in_val      <- candidate_idx[j]
                trial_score <-
                  compute_score(idx = c(subset_minus_out, in_val),
                                traits_mat = traits_mat,
                                div_fun =div_fun_internal,
                                metric = metric)

                # Track the best improvement found so far
                if (trial_score > best_overall_score) {
                  best_overall_score <- trial_score
                  best_out_val <- out_val
                  best_in_val <- in_val
                  best_out_pos_in_swappable <- i
                  best_in_pos_in_candidate  <- j
                }
              }
            }

            # Check if an improvement was actually found in this pass
            if (!is.null(best_out_val)) {
              # Update current collection
              current_idx[match(best_out_val, current_idx)] <- best_in_val
              current_score <- best_overall_score

              # Update indices in-place
              swappable_idx[best_out_pos_in_swappable] <- best_in_val
              candidate_idx[best_in_pos_in_candidate] <- best_out_val
              improved <- TRUE

              if (SampleCore.debug) {
                message(sprintf("Best-improvement | Iter %d: Swapped out %d for %d. New score: %f",
                                iter_1opt, best_out_val, best_in_val,
                                current_score))
              }
            }

          } else {

            ### First-improvement strategy ----

            # Nested loops
            for (i in seq_along(swappable_idx)) {
              out_val <- swappable_idx[i]
              subset_minus_out <- current_idx[current_idx != out_val]

              for (j in seq_along(candidate_idx)) {
                in_val <- candidate_idx[j]
                trial_score <- compute_score(c(subset_minus_out, in_val),
                                             traits_mat, div_fun_internal, metric)

                if (trial_score > current_score) {
                  # First improvement found - Apply swap immediately
                  current_idx[match(out_val, current_idx)] <- in_val
                  current_score <- trial_score

                  # swap in-place
                  swappable_idx[i] <- in_val # best_in enters swappable pool
                  candidate_idx[j] <- out_val # best_out enters candidate pool
                  improved <- TRUE

                  if (SampleCore.debug) {
                    message(sprintf("First-improvement | Iter %d: Swapped out %d for %d. New score: %f",
                                    iter_1opt, out_val, in_val, current_score))
                  }

                  break # Break inner loop
                }
              }
              if (improved) break # Break outer loop to restart 1-opt with new current_idx
            }
          }

          if (!improved) break # Local optimum reached
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


