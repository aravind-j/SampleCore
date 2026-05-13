#' Selection of Entries from Clusters/Groups by Random Sampling
#'
#' Select entries from cluster/groups in the entire collection by random
#' sampling according to allocation specified.
#'
#' For each cluster/group entries are selected randomly according to the
#' allocation provided
#' \insertCite{brown_core_1989,van_hintum_core_2000}{SampleCore}. Entries listed
#' as \code{always.selected} are mandatorily included in the selection. Warnings
#' are issued if requested allocation is smaller than the number of
#' always-selected entries in a cluster/group and/or when the cluster/group does
#' not contain enough remaining entries to fulfill the allocation.
#'
#' @template general-arg
#' @template sel-arg
#'
#' @returns A named list where each element contains the selected entry
#'   identifiers for a cluster/group.
#' @export
#'
#' @references
#'
#' \insertAllCited
#'
#' @examplesIf requireNamespace("cluster", quietly = TRUE)
#'
#' library(cluster)
#'
#' # Get data
#' data(cassava_EC_gp)
#'
#' data <- cbind(genotypes = rownames(cassava_EC_gp), cassava_EC_gp)
#' row.names(data) <- NULL
#'
#' # Prepare inputs
#' counts <- c(I = 31, II = 31, III = 18, IV = 35, V = 40, VI = 17)
#'
#' mand_accns <-
#'   c("TMe-34", "TMe-3423", "TMe-2018", "TMe-801", "TMe-551")
#'
#' # Specify the seed
#' set.seed(123)
#'
#' # Fetch selected accessions
#' sel_random_out <-
#'   select.random(data = data, names = "genotypes",
#'                 group = "Cluster", alloc = counts,
#'                 always.selected = mand_accns)
#'
#' sel_random_out
#'
#' # Get distance matrix - Only for visualization
#' quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW",
#'            "AVPW", "ARSR", "SRDM")
#' qual <- c("CUAL", "LNGS", "PTLC", "DSTA", "LFRT", "LBTEF", "CBTR", "NMLB",
#'           "ANGB", "CUAL9M", "LVC9M", "TNPR9M", "PL9M", "STRP", "STRC",
#'           "PSTR")
#'
#' # Convert qualitative data columns to factor
#' cassava_EC_gp[, qual] <- lapply(cassava_EC_gp[, qual], as.factor)
#'
#' # Standardise quantitative data column
#' cassava_EC_gp[, quant] <- lapply(cassava_EC_gp[, quant], function(x) {
#'   scale(x)[, 1]
#' })
#'
#' gp_vec <- setNames(as.character(data[, "Cluster"]), data[, "genotypes"])
#'
#' # Get the Gower's distance matrix
#' dist_matrix <- daisy(x = cassava_EC_gp[, c(qual, quant)],
#'                      metric = "gower")
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_random_out, use.names = FALSE))
#'
select.random <- function(data, names, group, alloc,
                          always.selected = NULL) {

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

      # Handle alloc smaller than fixed set
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

        # Handle insufficient remaining pool
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

        if (n_rem == 0) {

          sampled_accns <- fixed_accns

        } else {

        sampled_accns <- c(fixed_accns,
                           sample(rem_accns, n_rem, replace = FALSE))

        }

      }

      sampled_accns

    })

  names(out) <- names(alloc)

  return(out)

}


