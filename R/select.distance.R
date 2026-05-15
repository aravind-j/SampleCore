### This file is part of 'SampleCore' package for R.

### Copyright (C) 2024-2026, ICAR-NBPGR.
#
# SampleCore is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# SampleCore is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.r-project.org/Licenses/

#' Selection of Entries from Clusters/Groups on the basis of Genetic Distances
#'
#' Select entries from cluster/groups in the entire collection by genetic
#' distance based sampling according to allocation specified. \loadmathjax
#'
#' For each cluster/group, entries are selected by several methods from
#' within-cluster/group genetic distances between accessions according to the
#' allocation provided (See \strong{Methods}).
#'
#' Entries listed as \code{always.selected} are mandatorily included in the
#' selection. Warnings are issued if requested allocation is smaller than the
#' number of always-selected entries in a cluster/group and/or when the
#' cluster/group does not contain enough remaining entries to fulfill the
#' allocation.
#'
#' @section Methods:
#'
#' \subsection{Centrality Based Methods}{
#'
#'   Selects accessions that are most representative/closest to the
#'   cluster/group center.
#'
#'   \subsection{Medoid-like Representative Sampling by Minimal Mean
#'   Distance}{Selects medoid-like representatives as accessions with the
#'   smallest average distance to all others within the group
#'   \insertCite{kaufman_clustering_1987,kaufman_finding_1990}{SampleCore}.
#'
#'   For each accession \mjseqn{g}, the mean distance to all other accessions
#'   \mjseqn{h} is computed as:
#'
#'   \mjsdeqn{\bar{d}_g = \frac{1}{G} \sum_{h=1}^{G} d_{gh}}
#'
#'   Accessions are ranked by \mjseqn{\bar{d}_g} in ascending order and the top
#'   \mjseqn{n} are selected.
#'
#'   }
#'
#'   \subsection{Medoid-like Representative Sampling by Minimal Median
#'   Distance}{
#'
#'   Selects medoid-like representatives as accessions with the smallest median
#'   distance to all others within the group. This method is less influenced by
#'   outliers
#'   \insertCite{kaufman_clustering_1987,kaufman_finding_1990}{SampleCore}.
#'
#'   For each accession \mjseqn{g}, the median distance to all other accessions
#'   \mjseqn{h} is computed as:
#'
#'   \mjsdeqn{\tilde{d}_g = \text{median}_{h=1,\dots,G}(d_{gh})}
#'
#'   Accessions are ranked by \mjseqn{\tilde{d}_g} in ascending order and the
#'   top \mjseqn{n} are selected.
#'
#'   }
#'
#'   \subsection{Representative Sampling by Proximity to Group
#'   Centroid}{Selects accessions closest to the group centroid in principal
#'   coordinate space, computed via multivariate dispersion analysis using
#'   \code{\link[vegan]{betadisper}}
#'   \insertCite{anderson_distance-based_2006,anderson_multivariate_2006}{SampleCore}.
#'
#'   The distance of each accession \mjseqn{g} to the group centroid \mjseqn{C}
#'   in PCoA space is:
#'
#'   \mjsdeqn{\delta_g = \| \mathbf{p}_g - \mathbf{c} \|}
#'
#'   Where \mjseqn{\mathbf{p}_g} is the PCoA coordinate vector of accession
#'   \mjseqn{g} and \mjseqn{\mathbf{c}} is the group centroid. Accessions are
#'   ranked by \mjseqn{\delta_g} in ascending order and the top \mjseqn{n} are
#'   selected.
#'
#'   }
#'
#'   \subsection{Representative Sampling by Proximity to Group Spatial Median}{
#'   Selects accessions closest to the group spatial median in principal
#'   coordinate space, computed via multivariate dispersion analysis using
#'   \code{\link[vegan]{betadisper}}
#'   \insertCite{oneill_theory_2000}{SampleCore}.
#'
#'   The distance of each accession \mjseqn{g} to the group spatial median
#'   \mjseqn{M} is:
#'
#'   \mjsdeqn{\delta_g^* = \| \mathbf{p}_g - \mathbf{m} \|}
#'
#'   where \mjseqn{\mathbf{m}} is the spatial median of the group in PCoA
#'   space. Accessions are ranked by \mjseqn{\delta_g^*} in ascending order and
#'   the top \mjseqn{n} are selected.
#'
#'   }
#'
#' }
#'
#' \subsection{Peripheral/Extremity Based Methods}{
#'
#'   Selects accessions that are most dissimilar from the rest in a
#'   cluster/group i.e. the accessions which are in the boundary or outliers.
#'
#'   \subsection{Peripheral Sampling by Maximal Mean Distance}{Selects the most
#'   peripheral accessions as those with the largest average distance to all
#'   others within the group
#'   \insertCite{kaufman_clustering_1987,kaufman_finding_1990}{SampleCore}.
#'
#'   \mjsdeqn{\bar{d}_g = \frac{1}{G} \sum_{h=1}^{G} d_{gh}}
#'
#'   Accessions are ranked by \mjseqn{\bar{d}_g} in descending order and the
#'   top \mjseqn{n} are selected.
#'
#'   }
#'
#'   \subsection{Peripheral Sampling by Maximal Median Distance}{Selects the
#'   most peripheral accessions as those with the largest median distance to
#'   all others within the group
#'   \insertCite{kaufman_clustering_1987,kaufman_finding_1990}{SampleCore}.
#'
#'   \mjsdeqn{\tilde{d}_g = \text{median}_{h=1,\dots,G}(d_{gh})}
#'
#'   Accessions are ranked by \mjseqn{\tilde{d}_g} in descending order and the
#'   top \mjseqn{n} are selected.
#'
#'   }
#'
#'   \subsection{Peripheral Sampling by Maximal Eccentricity}{Selects
#'   accessions with the largest eccentricity — the maximum distance to any
#'   other accession in the group
#'   \insertCite{hage_eccentricity_1995}{SampleCore}.
#'
#'   \mjsdeqn{e_g = \max_{h=1,\dots,G} d_{gh}}
#'
#'   Accessions are ranked by \mjseqn{e_g} in descending order and the top
#'   \mjseqn{n} are selected. Eccentricity captures the worst-case
#'   dissimilarity of an accession rather than its average behaviour.
#'
#'   }
#'
#'   \subsection{Peripheral Sampling by Maximal Farness Centrality}{Selects
#'   accessions with the greatest total distance to all others, i.e. those most
#'   remote from the rest of the group
#'   \insertCite{sabidussi_centrality_1966}{SampleCore}.
#'
#'   \mjsdeqn{f_g = \sum_{h=1}^{G} d_{gh}}
#'
#'   Accessions are ranked by \mjseqn{f_g} in descending order and the top
#'   \mjseqn{n} are selected. Farness centrality is proportional to
#'   \mjseqn{\bar{d}_g} and differs from \code{mean.peripheral} only in that it
#'   uses the raw sum rather than the mean, producing identical rankings.
#'
#'   }
#'
#' }
#'
#' \subsection{Space-Filling/Coverage Methods}{
#'
#'   Select accessions that are spread maximally across the feature space in a
#'   cluster/group i.e. diversity sampling.
#'
#'   \subsection{Space-Filling Sampling via the Kennard-Stone
#'   Algorithm}{Selects \mjseqn{n} accessions that maximally and uniformly
#'   cover the distance space via the Kennard-Stone algorithm
#'   \insertCite{kennard_computer_1969}{SampleCore} (See
#'   \code{\link[prospectr]{kenStone}}).
#'
#'   Starting from the pair of accessions with the largest pairwise distance:
#'
#'   \mjsdeqn{\lbrace g_1, g_2 \rbrace = \underset{g,h}{\arg\max}\, d_{gh}}
#'
#'   each subsequent accession \mjseqn{g_k} is selected by maximising its
#'   minimum distance to the already-selected set \mjseqn{S}:
#'
#'   \mjsdeqn{g_k = \underset{g \notin S}{\arg\max} \min_{s \in S} d_{gs}}
#'
#'   This greedy procedure ensures even space coverage without relying on
#'   cluster structure.
#'
#'   }
#'
#'   \subsection{Space-Filling Sampling via the DUPLEX Algorithm}{Extends the
#'   Kennard-Stone algorithm to simultaneously construct a model set and a test
#'   set with similar distributions
#'   \insertCite{kennard_computer_1969,snee_validation_1977}{SampleCore}
#'   (\link[prospectr]{duplex}). Accessions are selected using Mahalanobis
#'   distance:
#'
#'   \mjsdeqn{d_M(g, h) = \sqrt{(\mathbf{x}_g - \mathbf{x}_h)^\top \Sigma^{-1}
#'   (\mathbf{x}_g - \mathbf{x}_h)}}
#'
#'   where \mjseqn{\Sigma} is the covariance matrix. At each step, the pair
#'   maximising \mjseqn{d_M} is split alternately between model and test sets,
#'   ensuring both sets span the full feature space.
#'
#'   }
#'
#'   \subsection{Space-Filling Sampling via the Honigs Algorithm}{Selects
#'   \mjseqn{n} accessions sequentially by maximising dissimilarity to the
#'   already-selected set \insertCite{honigs_unique-sample_1985}{SampleCore}
#'   (\link[prospectr]{honigs})
#'
#'   At each step \mjseqn{k}, the accession \mjseqn{g_k} maximising total
#'   distance to all previously selected accessions \mjseqn{S} is chosen:
#'
#'   \mjsdeqn{g_k = \underset{g \notin S}{\arg\max} \sum_{s \in S} d_{gs}}
#'
#'   This favours accessions that are collectively most dissimilar to the
#'   current selection, producing broad coverage of the distance space.
#'
#'   }.
#'
#'   \subsection{Space-Filling Sampling via Farthest-Point (Max-Min)
#'   Algorithm}{Selects \mjseqn{n} accessions by iteratively maximising the
#'   minimum distance to the current selected set — also known as the
#'   max-min or farthest-point sampling algorithm
#'   \insertCite{gonzalez_clustering_1985,dyer_simple_1985,hochbaum_best_1985}{SampleCore}.
#'
#'   \mjsdeqn{g_k = \underset{g \notin S}{\arg\max} \min_{s \in S} d_{gs}}
#'
#'   This is equivalent to Kennard-Stone but without the symmetric
#'   initialisation step. It provides a deterministic, greedy approximation to
#'   the \mjseqn{k}-centre problem:
#'
#'   \mjsdeqn{\min_{S \subset G,\, |S|=n} \max_{g \in G} \min_{s \in S} d_{gs}}
#'
#'   }
#'
#' }
#'
#' \subsection{Density Based Methods}{
#'
#'   Select points based on local neighbourhood density.
#'
#'   \subsection{Density-Based Sampling by Minimal Nearest-Neighbour
#'   Distance}{Selects accessions residing in the densest regions of the
#'   distance space, identified as those with the smallest nearest-neighbour
#'   distance
#'   \insertCite{cover_nearest_1967,fix_discriminatory_1989}{SampleCore}.
#'
#'   For each accession \mjseqn{g}, the nearest-neighbour distance is:
#'
#'   \mjsdeqn{\text{nn}_g = \min_{h \neq g} d_{gh}}
#'
#'   Accessions are ranked by \mjseqn{\text{nn}_g} in ascending order and the
#'   top \mjseqn{n} are selected. Small \mjseqn{\text{nn}_g} indicates that
#'   \mjseqn{g} resides in a dense cluster; this method preferentially samples
#'   from high-density regions.
#'
#'   }
#'
#' }
#'
#' \subsection{Cluster Based Methods}{
#'
#'   These methods partition the cluster/group space into sub-clusters/groups,
#'   then samples from each one.
#'
#'   \subsection{Globally Optimal Medoid Sampling via Partitioning Around
#'   Medoids (PAM)}{
#'
#'   Selects a set of \mjseqn{n} medoids that jointly minimise the total
#'   distance of every accession to its nearest medoid, via
#'   \code{\link[cluster]{pam}}.
#'
#'   The objective function minimised is:
#'
#'   \mjsdeqn{\min_{S \subset G,\, |S|=n} \sum_{g=1}^{G} \min_{s \in S} d_{gs}}
#'
#'   Unlike \code{"mean.medoid"}, medoids are co-optimised as a set, ensuring
#'   they collectively represent the full distribution of the group rather than
#'   independently scoring each accession.
#'
#'   }
#'
#'   \subsection{Cluster-Based Sampling via K-means (Naes Method)}{Partitions
#'   accessions into \mjseqn{n} clusters via k-means applied to the distance
#'   matrix (See \code{\link[prospectr]{naes}}), then selects the accession
#'   closest to each cluster centre as the representative
#'   \insertCite{naes_design_1987,naes_user-friendly_2017}{SampleCore}.
#'
#'   The k-means objective minimised is:
#'
#'   \mjsdeqn{\min \sum_{k=1}^{n} \sum_{g \in C_k} d_{g, \mu_k}^2}
#'
#'   where \mjseqn{C_k} is the \mjseqn{k}-th cluster and \mjseqn{\mu_k} is its
#'   centre. One representative per cluster is returned, ensuring broad,
#'   partition-aware coverage.
#'
#'   }
#'
#'   \subsection{Cluster-Based Sampling via Hierarchical Clustering with
#'   Random Selection}{Partitions accessions into \mjseqn{n} clusters by
#'   cutting a hierarchical clustering dendrogram at height \mjseqn{k = n},
#'   then randomly samples one accession from each cluster
#'   \insertCite{ward_Hierarchical_1963,li_studies_2002}{SampleCore}.
#'
#'   The dendrogram is built by agglomerative hierarchical clustering using the
#'   linkage criterion specified by \code{\link[stats]{hclust}}. For
#'   clusters \mjseqn{C_1, \dots, C_n}, one accession is drawn uniformly at
#'   random from each:
#'
#'   \mjsdeqn{g_k \sim \text{Uniform}(C_k), \quad k = 1, \dots, n}
#'
#'   This introduces stochasticity within a structured partition, balancing
#'   coverage with randomness.
#'
#'   }
#'
#'   \subsection{Cluster-Based Sampling via Hierarchical Clustering with Medoid
#'   Selection}{Partitions accessions into \mjseqn{n} clusters by cutting a
#'   hierarchical clustering dendrogram at height \mjseqn{k = n}, then selects
#'   the within-cluster medoid as the representative of each cluster
#'   \insertCite{kaufman_clustering_1987,ward_Hierarchical_1963}{SampleCore}.
#'
#'   For each cluster \mjseqn{C_k}, the medoid is the accession minimising
#'   total within-cluster distance:
#'
#'   \mjsdeqn{g_k^* = \underset{g \in C_k}{\arg\min} \sum_{h \in C_k} d_{gh}}
#'
#'   This combines the structured partitioning of hierarchical clustering with
#'   deterministic, centrality-based representative selection.
#'
#'   }
#'
#' }
#'
#' @template general-arg
#' @template sel-arg
#' @template dist-arg
#' @template seldist-arg
#' @param method The method for sampling accessions from each cluster/group.
#'   Either \code{"mean.medoid"}, \code{"median.medoid"},
#'   \code{"nearest.centroid"}, \code{"nearest.median"},
#'   \code{"mean.peripheral"}, \code{"median.peripheral"},
#'   \code{"eccentricity"}, \code{"farness.centrality"}, \code{"kennard.stone"},
#'   \code{"duplex"}, \code{"honigs"}, \code{"farthest.point"},
#'   \code{"nearest.neighbour"}, \code{"naes"}, \code{"optim.medoid"},
#'   \code{"hclust.random"} or \code{"hclust.medoid"}. See \strong{Methods}.
#'
#' @returns A named list where each element contains the selected entry
#'   identifiers for a cluster/group.
#'
#' @importFrom cluster pam
#' @importFrom graphics par title
#' @importFrom grDevices n2mfrow
#' @importFrom MASS isoMDS
#' @importFrom prospectr duplex honigs kenStone naes
#' @importFrom stats cmdscale cutree hclust
#' @importFrom vegan betadisper
#'
#' @export
#'
#' @references
#'
#' \insertAllCited
#'
#' @seealso \code{\link[SampleCore]{select.random}},
#'   \code{\link[SampleCore]{select.diversity}}
#'
#' @examplesIf requireNamespace("cluster", quietly = TRUE) & requireNamespace("ggplot2", quietly = TRUE)
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Prepare example data
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#'
#' library(cluster)
#' library(ggplot2)
#'
#' data(cassava_EC_gp)
#'
#' set.seed(123)
#' data <- cassava_EC_gp[sample(1:nrow(cassava_EC_gp), 500), ]
#'
#' quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW", "AVPW",
#'            "ARSR", "SRDM")
#' qual <- c("CUAL", "LNGS", "PTLC", "DSTA", "LFRT", "LBTEF", "CBTR", "NMLB",
#'           "ANGB", "CUAL9M", "LVC9M", "TNPR9M", "PL9M", "STRP", "STRC",
#'           "PSTR")
#'
#' data[, qual] <- lapply(data[, qual], as.factor)
#'
#' # Get the Gower's distance matrix
#' dist_matrix <- daisy(x = data[, c(qual, quant)],
#'                      metric = "gower")
#'
#'
#' data <- cbind(genotypes = rownames(data), data)
#' row.names(data) <- NULL
#'
#' # Prepare inputs
#' counts <- c(I = 16, II = 15, III = 9, IV = 18, V = 20, VI = 8)
#'
#' mand_accns <-
#'   c("TMe-2018", "TMe-801", "TMe-3191", "TMe-1830", "TMe-1790")
#'
#' gp_vec <- setNames(as.character(data[, "Cluster"]), data[, "genotypes"])
#'
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Fetch selected accessions by centrality based methods
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' # Medoid-like Representative Sampling by Minimal Mean Distance
#' sel_mean_medoid_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "mean.medoid")
#' sel_mean_medoid_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_mean_medoid_out,
#'           use.names = FALSE)) +
#'   labs(title = "mean.medoid")
#'
#' # Medoid-like Representative Sampling by Minimal Median Distance
#' sel_median_medoid_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "median.medoid")
#' sel_median_medoid_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_median_medoid_out,
#'           use.names = FALSE)) +
#'   labs(title = "median.medoid")
#'
#' # Representative Sampling by Proximity to Group Centroid
#' sel_group_centroid_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "nearest.centroid")
#' sel_group_centroid_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_group_centroid_out,
#'           use.names = FALSE)) +
#'   labs(title = "nearest.centroid")
#'
#' # Representative Sampling by Proximity to Group Spatial Median
#' sel_group_median_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "nearest.median")
#' sel_group_median_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_group_median_out,
#'           use.names = FALSE)) +
#'   labs(title = "nearest.median")
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Fetch selected accessions by peripheral/extremity based methods
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' # Peripheral Sampling by Maximal Mean Distance
#' sel_mean_peripheral_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "mean.peripheral")
#' sel_mean_peripheral_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_mean_peripheral_out,
#'           use.names = FALSE)) +
#'   labs(title = "mean.peripheral")
#'
#' # Peripheral Sampling by Maximal Median Distance
#' sel_median_peripheral_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "median.peripheral")
#' sel_median_peripheral_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_median_peripheral_out,
#'           use.names = FALSE)) +
#'   labs(title = "median.peripheral")
#'
#' # Peripheral Sampling by Maximal Eccentricity
#' sel_eccentricity_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "eccentricity")
#' sel_eccentricity_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_eccentricity_out,
#'           use.names = FALSE)) +
#'   labs(title = "eccentricity")
#'
#' # Peripheral Sampling by Maximal Farness Centrality
#' sel_far_cent_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "farness.centrality")
#' sel_far_cent_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_far_cent_out,
#'           use.names = FALSE)) +
#'   labs(title = "farness.centrality")
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Fetch selected accessions by space-Filling/coverage methods
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' # Space-Filling Sampling via the Kennard-Stone Algorithm
#' sel_ks_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "kennard.stone")
#' sel_ks_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_ks_out,
#'           use.names = FALSE)) +
#'   labs(title = "kennard.stone")
#'
#' # Space-Filling Sampling via the DUPLEX Algorithm
#' sel_duplex_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "duplex")
#' sel_duplex_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_duplex_out,
#'           use.names = FALSE)) +
#'   labs(title = "duplex")
#'
#' # Space-Filling Sampling via the Honigs Algorithm
#' sel_honigs_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "honigs")
#' sel_honigs_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_honigs_out,
#'           use.names = FALSE)) +
#'   labs(title = "honigs")
#'
#' # Space-Filling Sampling via Farthest-Point (Max-Min) Algorithm
#' sel_far_pt_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "farthest.point")
#' sel_far_pt_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_far_pt_out,
#'           use.names = FALSE)) +
#'   labs(title = "farthest.point")
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Fetch selected accessions by density based methods
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' # Density-Based Sampling by Minimal Nearest-Neighbour Distance
#' sel_nn_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "nearest.neighbour")
#' sel_nn_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_nn_out,
#'           use.names = FALSE)) +
#'   labs(title = "nearest.neighbour")
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Fetch selected accessions by cluster based methods
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' # Globally Optimal Medoid Sampling via Partitioning Around Medoids (PAM)
#' sel_pam_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "optim.medoid")
#' sel_pam_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_pam_out,
#'           use.names = FALSE)) +
#'   labs(title = "optim.medoid")
#'
#' # Cluster-Based Sampling via K-means (Naes Method)
#' sel_naes_out <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "naes")
#' sel_naes_out
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_naes_out,
#'           use.names = FALSE)) +
#'   labs(title = "naes")
#'
#' # Cluster-Based Sampling via Hierarchical Clustering with Random Selection
#'
#' ## UPGMA
#' sel_hclust_random_out1 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.random",
#'                   hclust.method = "average")
#' sel_hclust_random_out1
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_random_out1,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.random", subtitle = "average")
#'
#' ## Single-linkage
#' sel_hclust_random_out2 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.random",
#'                   hclust.method = "single")
#' sel_hclust_random_out2
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_random_out2,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.random", subtitle = "single")
#'
#' ## Complete-linkage
#' sel_hclust_random_out3 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.random",
#'                   hclust.method = "complete")
#' sel_hclust_random_out3
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_random_out3,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.random", subtitle = "complete")
#'
#' ## Ward's D
#' sel_hclust_random_out4 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.random",
#'                   hclust.method = "ward.D")
#' sel_hclust_random_out4
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_random_out4,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.random", subtitle = "ward.D")
#'
#' ## WPGMA
#' sel_hclust_random_out5 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.random",
#'                   hclust.method = "mcquitty")
#' sel_hclust_random_out5
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_random_out5,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.random", subtitle = "mcquitty")
#'
#' ## WPGMC
#' sel_hclust_random_out6 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.random",
#'                   hclust.method = "median")
#' sel_hclust_random_out6
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_random_out6,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.random", subtitle = "median")
#'
#' ## UPGMC
#' sel_hclust_random_out7 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.random",
#'                   hclust.method = "centroid")
#' sel_hclust_random_out7
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_random_out7,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.random", subtitle = "centroid")
#'
#' ## Ward's D2
#' sel_hclust_random_out8 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.random",
#'                   hclust.method = "ward.D2")
#' sel_hclust_random_out8
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_random_out8,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.random", subtitle = "ward.D2")
#'
#' # Cluster-Based Sampling via Hierarchical Clustering with Medoid Selection
#'
#' ## UPGMA
#' sel_hclust_medoid_out1 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.medoid",
#'                   hclust.method = "average")
#' sel_hclust_medoid_out1
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_medoid_out1,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.medoid", subtitle = "average")
#'
#' ## Single-linkage
#' sel_hclust_medoid_out2 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.medoid",
#'                   hclust.method = "single")
#' sel_hclust_medoid_out2
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_medoid_out2,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.medoid", subtitle = "single")
#'
#' ## Complete-linkage
#' sel_hclust_medoid_out3 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.medoid",
#'                   hclust.method = "complete")
#' sel_hclust_medoid_out3
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_medoid_out3,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.medoid", subtitle = "complete")
#'
#' ## Ward's D
#' sel_hclust_medoid_out4 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.medoid",
#'                   hclust.method = "ward.D")
#' sel_hclust_medoid_out4
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_medoid_out4,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.medoid", subtitle = "ward.D")
#'
#' ## WPGMA
#' sel_hclust_medoid_out5 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.medoid",
#'                   hclust.method = "mcquitty")
#' sel_hclust_medoid_out5
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_medoid_out5,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.medoid", subtitle = "mcquitty")
#'
#' ## WPGMC
#' sel_hclust_medoid_out6 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.medoid",
#'                   hclust.method = "median")
#' sel_hclust_medoid_out6
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_medoid_out6,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.medoid", subtitle = "median")
#'
#' ## UPGMC
#' sel_hclust_medoid_out7 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.medoid",
#'                   hclust.method = "centroid")
#' sel_hclust_medoid_out7
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_medoid_out7,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.medoid", subtitle = "centroid")
#'
#' ## Ward's D2
#' sel_hclust_medoid_out8 <-
#'   select.distance(data = data, names = "genotypes",
#'                   group = "Cluster", alloc = counts,
#'                   dist.mat = dist_matrix,
#'                   always.selected = mand_accns,
#'                   method = "hclust.medoid",
#'                   hclust.method = "ward.D2")
#' sel_hclust_medoid_out8
#'
#' plot_dist(d = dist_matrix, method = "isomds",
#'           gp = gp_vec,
#'           highlight =  unlist(sel_hclust_medoid_out8,
#'           use.names = FALSE)) +
#'   labs(title = "hclust.medoid", subtitle = "ward.D2")
#'
select.distance <- function(data, names, group, alloc,
                            dist.mat, always.selected = NULL,
                            method = c("mean.medoid", "median.medoid",
                                       "nearest.centroid", "nearest.median",
                                       "mean.peripheral", "median.peripheral",
                                       "eccentricity", "farness.centrality",
                                       "kennard.stone", "duplex", "honigs",
                                       "farthest.point",
                                       "nearest.neighbour",
                                       "naes", "optim.medoid",
                                       "hclust.random", "hclust.medoid"),
                            hclust.method = c("average", "single",
                                              "complete",
                                              "ward.D", "mcquitty",
                                              "median", "centroid",
                                              "ward.D2")) {

  SampleCore.debug <- getOption("SampleCore.debug", default = FALSE)

  method <- match.arg(method)

  checks.sample.core(data = data, size = NULL,
                     names = names, group = group,
                     dist.mat = dist.mat,
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

  if (SampleCore.debug) {
    par(mfrow = n2mfrow(length(alloc)))
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

        avail_rem <- length(rem_accns)

        if (avail_rem < n_rem) {

          warning(
            sprintf(
              paste0('Group %s has only %d additional accessions but ',
                     '"alloc" requires %d additional entries. ',
                     'Taking all available accessions.'),
              g, avail_rem, n_rem
            ),
            call. = FALSE
          )

          n_rem <- avail_rem
        }

        if (n_rem == 1) {
          sampled_accns <- rem_accns[which.min(avg_dist)]
        } else if (avail_rem < 3 && method %in% c("kennard.stone", "duplex",
                                                  "honigs", "naes")) {
          warning(
            sprintf(
              paste0('Group %s: method "%s" requires at least 3 accessions ',
                     'but only %d available. Falling back to "mean.medoid".'),
              g, method, avail_rem
            ),
            call. = FALSE
          )
          center_idx <- order(avg_dist)[1:n_rem]
          sampled_accns <- rem_accns[center_idx]

        } else { # Method dispatch

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

          ## Max-min / Farthest-point sampling ----
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
            tree_out <- hclust(as.dist(sub_d), method = hclust.method)
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
        }
      }

      return(c(sampled_accns, fixed_accns))

    })

  names(out) <- names(alloc)

  if (SampleCore.debug) {

    gp_vec <- setNames(as.character(data[, group]), data[, names])

    viz_out <-
      plot_dist(d = dist.mat, method = "isomds", gp = gp_vec,
              highlight = unlist(out, use.names = FALSE))
    print(viz_out)

  }

  return(out)
}

## Max-min / Farthest-point sampling function ----
farthest_sampling <- function(d, n_select,
                              start = sample(nrow(d), 1)) {

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
