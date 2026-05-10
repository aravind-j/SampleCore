#' Selection of Entries from Clusters/Groups by Random Sampling
#'
#' Select entries from cluster/groups in the entire collection by random
#' sampling according to allocation specified.
#'
#' For each cluster/group entries are selected randomly according to the
#' allocation provided \insertCite{brown_core_1989;
#' van_hintum_core_2000}{SampleCore}. Entries listed as \code{always.selected}
#' are mandatorily included in the selection. Warnings are issued if requested
#' allocation is smaller than the number of always-selected entries in a
#' cluster/group and/or when the cluster/group does not contain enough remaining
#' entries to fulfill the allocation.
#'
#' @template general-arg
#' @template sel-arg
#'
#' @returns A named list where each element contains the selected entry
#'   identifiers for a cluster/group.
#' @export
#'
#' @examples
select.random <- function(data, names, group, alloc, always.selected) {

  checks.sample.core(data = data, size = NULL,
                     names = names, group = group,
                     dist.mat = NULL,
                     quantitative = NULL,
                     qualitative = NULL,
                     log.base = NULL,
                     alloc = alloc,
                     always.selected = always.selected,
                     mode = "sel")

  out <-
    lapply(names(alloc), function(g) {

      group_accns <- data[data[[group]] == g, names]

      # always-selected values belonging to this group
      fixed_accns <- intersect(always.selected, group_accns)

      requested_n <- alloc[[g]]

      # Remaining pool after fixed selections
      rem_accns <- setdiff(group_accns, fixed_accns)

      n_rem <- requested_n - length(fixed_accns)

      if (n_rem < 0) {

        warning(
          sprintf(
            paste0('Group %s: "alloc" (%d) is smaller than number of ',
                   '"always.selected" values (%d). ',
                   'Taking only "always.selected" values.'),
            g, requested_n, length(fixed_accns)
          ),
          call. = FALSE
        )

        sampled_accns <- fixed_accns
      } else {

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

        sampled_accns <- c(fixed_accns, sample(rem_accns, n_rem))

      }

      sampled_accns

    })

  names(out) <- names(alloc)

  return(out)

}


