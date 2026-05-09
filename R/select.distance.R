#' Selection of Entries from Clusters/Groups on the basis of Genetic Distances
#'
#' Select entries from cluster/groups in the entire collection by genetic
#' distance based sampling according to allocation specified.
#'
#' For each cluster/group, entries are selected on the basis of several metrics
#' estimated from the within cluster/group genetic distances between accessions
#' according to the allocation provided.
#'
#' \subsection{Centrality Based Methods}{
#'
#'   \subsection{Mean Medoid}{medoid-like representatives of the dataset as accessions with smallest average distance to all others}
#'
#'   \subsection{Median Medoid}{.... smallest average distance.... Less influenced by outliers }
#'
#'   \subsection{Nearest to centroid}{}
#'
#'   \subsection{Nearest to median}{}
#'
#' }
#'
#' \subsection{Peripheral/Extremity Based Methods}{
#'
#'   \subsection{Mean Peripheral}{}
#'
#'   \subsection{Median Peripheral}{}
#'
#'   \subsection{Eccentricity}{}
#'
#'   \subsection{Farness Centrality}{}
#'
#' }
#'
#' \subsection{Space-Filling/Coverage Methods}{
#'
#'   \subsection{Kennard-Stone Algorithm}{}
#'
#'   \subsection{DUPLEX Algorithm}{}
#'
#'   \subsection{Farness centrality}{}
#'
#'   \subsection{Honigs Algorithm}{}
#'
#'   \subsection{Max-Min/Farthest-Point Sampling}{}
#'
#' }
#'
#' \subsection{Density Based Methods}{}
#'
#' \subsection{Cluster Based Methods}{
#'
#'   \subsection{PAM}{}
#'
#' }
#'
#' Entries listed as \code{always.selected} are mandatorily included in the
#' selection. Warnings are issued if requested allocation is smaller than the
#' number of always-selected entries in a cluster/group and/or when the
#' cluster/group does not contain enough remaining entries to fulfill the
#' allocation.
#'
#' @template general-arg
#' @template sel-arg
#' @template seldist-arg
#'
#' @returns A named list where each element contains the selected entry
#'   identifiers for a cluster/group.
#' @export
#'
#' @examples
select.distance <- function(data, names, group, alloc,
                            dist.mat, always.selected,
                            method = c("mean.medoid", "median.medoid",
                                       "nearest.centroid", "nearest.median",
                                       "mean.peripheral", "median.peripheral",
                                       "eccentricity", "farness.centrality",
                                       "kennard.stone", "duplex", "honigs",
                                       "farthest.point",
                                       "nearest.neighbour",
                                       "naes", "optim.medoid",
                                       "hclust.random", "hclust.medoid",
                                       "inv.tocher"),
                            hclust.method = c("average", "single",
                                              "complete",
                                              "ward.D", "mcquitty",
                                              "median", "centroid",
                                              "ward.D2")) {

  method <- match.arg(method)

  checks.sample.core(data = data, size = NULL,
                     names = names, group = group,
                     dist.mat = dist.mat,
                     quantitative = NULL,
                     qualitative = NULL,
                     log.base = NULL,
                     alloc = alloc,
                     always.selected = always.selected,
                     mode = "sel")

  # Prepare data ----

  dmat <- as.matrix(dist.mat)

  gp_memb <- data[, group]
  dist_labels <- labels(dist.mat)

  # Align labels
  names(gp_memb) <- dist_labels
  gp_memb <- gp_memb[rownames(dmat)]

  if (method == "nearest.centroid") {
    bdout <- vegan::betadisper(d = dist.mat, group = data[, group],
                               type = "centroid", sqrt.dist = TRUE)
    dist_to_centroid <- split(bdout$distances, bdout$group)
  }

  if (method == "nearest.median") {
    bdout <- vegan::betadisper(d = dist.mat, group = data[, group],
                               type = "median", sqrt.dist = TRUE)
    dist_to_median <- split(bdout$distances, bdout$group)
  }

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


        idx <- which(labels(dist.mat) %in% rem_accns)

        sub_d <- dmat[idx, idx]

        avg_dist <- rowMeans(as.matrix(sub_d))
        med_dist <- apply(as.matrix(sub_d), 1, median)

        if (method %in% c("mean.centroid", "median.centroid")) {
          is_group <- bdout$group == g
          group_distances <- bdout$distances[is_group]

          group_distances[names(group_distances) %in% rem_accns]
        }

        avail_rem <- length(rem_accns)

        if (avail_rem < n_rem) {

          warning(
            sprintf(
              'Group %s has only %d additional accessions but "alloc" ',
              'requires %d additional entries. ',
              'Taking all available accessions.',
              g, avail_rem, n_rem
            ),
            call. = FALSE
          )

          n_rem <- avail_rem
        }

        # Centrality Based ----

        ## Mean medoid ----
        if (method == "mean.medoid") {
          center_idx <- order(avg_dist)[1:n_rem]
          sampled_accns <- rem_accns[center_idx]
        }

        ## Median medoid ----
        if (method == "median.medoid") {
          center_idx2 <- order(med_dist)[1:n_rem]
          sampled_accns <- rem_accns[center_idx2]
        }

        ## Nearest to centroid ----
        if (method == "nearest.centroid") {
          dist_to_centroid_rem <- dist_to_centroid[[g]][rem_accns]
          center_idx3 <- order(dist_to_centroid_rem)[1:n_rem]
          sampled_accns <- names(dist_to_centroid[[g]])[center_idx3]
        }

        ## Nearest to centroid ----
        if (method == "nearest.median") {
          dist_to_median_rem <- dist_to_median[[g]][rem_accns]
          center_idx4 <- order(dist_to_median_rem)[1:n_rem]
          sampled_accns <- names(dist_to_median[[g]])[center_idx4]
        }

        # Peripheral / Extremity-Based ----

        ## Mean peripheral -----
        if (method == "mean.peripheral") {
          extreme_idx <- order(avg_dist, decreasing = TRUE)[1:n_rem]
          sampled_accns <- rem_accns[extreme_idx]
        }

        ## Median peripheral ----
        if (method == "median.peripheral") {
          extreme_idx2 <- order(med_dist, decreasing = TRUE)[1:n_rem]
          sampled_accns <- rem_accns[extreme_idx2]
        }

        ## Eccentricity ----
        if (method == "eccentricity") {
          ecc <- apply(as.matrix(sub_d), 1, max)
          ecc_idx <- order(ecc, decreasing = TRUE)[1:n_rem]
          sampled_accns <- rem_accns[ecc_idx]
        }

        ## Farness Centrality ----
        if (method == "farness.centrality") {
          far <- apply(as.matrix(sub_d), 1, sum)
          far_idx <- order(far, decreasing = TRUE)[1:n_rem]
          sampled_accns <- rem_accns[far_idx]
        }

        # Space-Filling / Coverage Methods ----

        ## Kennard-Stone algorithm ----
        if (method == "kennard.stone") {
          kens_res <- prospectr::kenStone(sub_d, k = n_rem)
          sampled_accns <- rem_accns[kens_res$model]
        }

        ## DUPLEX algorithm ----
        if (method == "duplex") {
          dupl_res <- prospectr::duplex(sub_d, k = n_rem, metric = "mahal")
          sampled_accns <- rem_accns[dupl_res$model]
        }

        ## Honigs algorithm ----
        if (method == "honigs") {
          ho_res <- prospectr::honigs(sub_d, k = n_rem)
          sampled_accns <- rem_accns[ho_res$model]
        }

        ## Max-min / farthest-point sampling ----
        if (method == "farthest.point") {
          farpoint_idx <- farthest_sampling(d = sub_d, n_select = n_rem)
          sampled_accns <- rem_accns[farpoint_idx]
        }

        # Density Based Methods ----

        ## Nearest-neighbour (NN) ----
        if (method == "nearest.neighbour") {
          nn <- apply(as.matrix(sub_d) + diag(Inf, nrow(sub_d)), 1, min)
          nn_idx <- order(nn)[1:n_rem]
          sampled_accns <- rem_accns[nn_idx]
        }

        # Cluster-Based Methods ----

        ## Partitioning around medoids ---
        if (method == "optim.medoid") {
          pam_res <- cluster::pam(sub_d, k = n_rem, diss = TRUE)
          sampled_accns <- pam_res$medoids
        }

        ## K-means sampling (Naes Method) ----
        if (method == "naes") {
          kms_res <- prospectr::naes(sub_d, k = n_rem)
          sampled_accns <- rem_accns[kms_res$model]
        }

        ## Hierarchical Clustering with random selection ----
        if (method %in% c("hclust.random", "hclust.medoid")) {
          tree_out <- hclust(as.dist(sub_d), method = "average")
          k_out <- cutree(tree_out, k = n_rem)
          k_out_list <- split(names(k_out), k_out)
        }

        if (method == "hclust.random") {
          sampled_accns <- sapply(k_out_list,
                                  function(x) {
                                    sample(x, 1)
                                  })
          sampled_accns <- unname(sampled_accns)
        }

        ## Hierarchical Clustering and medoids ----
        if (method == "hclust.medoid") {
          sampled_accns <- sapply(k_out_list, function(members) {
            # cluster submatrix
            k_d <- dmat[members, members, drop = FALSE]
            # total within-cluster distance
            tot_dist <- rowSums(k_d)
            # medoid
            names(which.min(tot_dist))
          })
          sampled_accns <- unname(sampled_accns)
        }

        ## Inverse Tocher ----
        if (method == "inv.tocher") {
          invtoch_idx <- inverse_tocher(sub_d, n_rem)
          sampled_accns <- rem_accns[invtoch_idx]
        }

      }

      if (g == gp_memb[1]) {
        plot_dist(d = sub_d, method = "isomds", highlight = sampled_accns)
      }

      return(c(sampled_accns, fixed_accns))

    })

  names(out) <- names(alloc)

  return(out)
}


farthest_sampling <- function(d, n_select, start = 1) {

  d <- as.matrix(d)

  selected <- start

  while(length(selected) < n_select) {

    remaining <- setdiff(seq_len(nrow(d)), selected)

    min.dist <- sapply(remaining, function(i) {
      min(d[i, selected])
    })

    next.idx <- remaining[which.max(min.dist)]

    selected <- c(selected, next.idx)
  }

  selected
}


inverse_tocher <- function(d, n_select) {

  d <- as.matrix(d)
  n <- nrow(d)

  # most distant pair
  idx <- which(d == max(d), arr.ind = TRUE)[1, ]
  selected <- unique(idx)

  while(length(selected) < n_select) {

    remaining <- setdiff(1:n, selected)

    avgd <- sapply(remaining, function(i) {
      mean(d[i, selected])
    })

    next_idx <- remaining[which.max(avgd)]

    selected <- c(selected, next_idx)
  }

  selected
}


plot_dist <- function(d, method = c("cmds", "isomds", "tsne"),
                      highlight) {

  # Classical MDS
  if (method == "cmds") {
    fit <- cmdscale(d, k = 2)
  }

  # Non-metric MDS
  # when distances are non-Euclidean or rank information matters
  if (method == "isomds") {
    fit <- MASS::isoMDS(as.dist(d))
    fit <- fit$points
  }

  # t-SNE
  if (method == "tsne") {
    fit <- Rtsne::Rtsne(as.matrix(d), is_distance = TRUE)
    fit <- fit$Y
    rownames(fit) <- labels(d)[[1]]
  }

  fit <- data.frame(fit)

  cols <- rep("black", nrow(fit))
  cols[which(row.names(fit) %in% highlight)] <- "red"

  plot(fit[,1], fit[,2],
       # xlab = "MDS1",
       # ylab = "MDS2",
       xlab = "X",
       ylab = "Y",
       col = cols,
       pch = 19)
}

