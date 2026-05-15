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

#' Allocation of Entries to be Selected from Clusters/Groups based on
#' Distance-based Diversity Metrics for Core Collection Development
#'
#' Estimate the number of entries to be allocated from each cluster/group in
#' the entire collection to construct a core collection on the basis of
#' different metrics computed from within cluster/group distances. The following
#' strategies are implemented. \loadmathjax
#' \itemize{
#' \item{Diversity (Distance based)}
#' \item{Diversity (Distance based) & Proportional}
#' \item{Diversity (Distance based) & Logarithmic}
#' \item{Diversity (Distance based) & Square root}}
#'
#' @section Details:
#'
#'   The number of entries to be chosen from each cluster is estimated either on
#'   the basis of diversity of entries within that cluster/group alone or in
#'   combination with the size of the cluster/group (See \strong{Methods}).
#'
#'   The within-cluster/group diversity is estimated as several metrics from the
#'   within cluster/group genetic distances between accessions (See
#'   \strong{Metrics}).
#'
#'   \insertCite{franco_sampling_2005;textual}{SampleCore} proposed a method
#'   based on mean Gower's distance \insertCite{gower_general_1971}{SampleCore}
#'   which was also extended to other distance measure averages named D
#'   Allocation strategy \insertCite{franco_sampling_2006}{SampleCore}. These
#'   methods were also combined with the proportional and logarithmic methods.
#'   For example, the GP and GL strategy of
#'   \insertCite{bisht_assessment_1999;textual}{SampleCore} and
#'   \insertCite{mahajan_sampling_1999;textual}{SampleCore} as well as the NY
#'   and LD allocation methods of
#'   \insertCite{franco_sampling_2005;textual}{SampleCore}.
#'
#' @template divmethods-section
#'
#' @section Metrics:
#'
#'   \subsection{Summary/Decriptive statistics}{These include mean, median,
#'   maximum and range of genetic distances between entries in a cluster.}
#'
#'   \subsection{Mean nearest-neighbour distance (\mjseqn{MNND})}{It is the
#'   average, across all entries, of the distance to each entry’s
#'   closest other entry (\mjseqn{d_{g_{min}}}), based on a genetic given
#'   distance matrix \insertCite{clark_distance_1954}{SampleCore}.
#'
#'   For each entry, the nearest-neighbour distance (\mjseqn{d_{g_{min}}})
#'   is the smallest non-zero distance with any other entry.
#'
#'   \mjsdeqn{d_{g_{min}} = \min_{h \ne g} d_{gh}}
#'
#'   The Mean nearest-neighbour distance (\mjseqn{MNND}) can then be computed
#'   as:
#'
#'   \mjsdeqn{\textrm{MNND} = \frac{1}{G} \sum_{g=1}^{G} d_g}
#'
#'   Where, (\mjseqn{g}) is the index of an entry in a genetic distance
#'   matrix, \mjseqn{h} is the index of all other genotypes and \mjseqn{G} is
#'   the total number of genotypes in a cluster/group.
#'
#'   }
#'
#'   \subsection{Minimum spanning tree length (\mjseqn{MSTL})}{It is defined as
#'   the sum of edge weights in the minimum spanning tree constructed from the
#'   genetic distance matrix of entries within a cluster/group.  A minimum
#'   spanning tree (MST) connects all entries such that the total distance
#'   is minimized and no cycles are formed. It represents the most efficient
#'   way to connect all entries based on pairwise genetic distances
#'   \insertCite{gower_minimum_1969}{SampleCore}.
#'
#'   For genetic distance \mjseqn{d_{gh}} between entries \mjseqn{g} and
#'   \mjseqn{h}, the MST is a subset of edges that connects all
#'   \mjseqn{G} entries with exactly \mjseqn{G - 1} edges and minimum total
#'   weight. The MST length (\mjseqn{MSTL}) can then be computed as:
#'
#'   \mjsdeqn{\textrm{MSTL} = \sum_{(g,h) \in \mathcal{T}} d_{gh}}
#'
#'   Where \mjseqn{\mathcal{T}} denotes the set of edges in the MST.
#'
#'   }
#'
#'   \subsection{Mean distance to centroid and median (\mjseqn{MDC},
#'   \mjseqn{MDM})}{These quantify the average dispersion of entries within a
#'   cluster/group relative to a central point in multivariate space derived
#'   from the genetic distance matrix.
#'
#'   The centroid represents the multivariate mean position of all entries
#'   in a cluster
#'   \insertCite{sokal_principles_1963,sneath_numerical_1973}{SampleCore}.,
#'   whereas the median (spatial median) provides a robust central location
#'   that is less influenced by extreme values
#'   \insertCite{bradley_constrained_1999}{SampleCore}.
#'
#'   For \mjseqn{d_{gC}} and \mjseqn{d_{gM}} distances of entry
#'   \mjseqn{g} from the centroid \mjseqn{C} and median \mjseqn{M},
#'   respectively. These measures are computed as:
#'
#'   \mjsdeqn{\textrm{MDC} = \frac{1}{G} \sum_{g=1}^{G} d_{gC}}
#'
#'   \mjsdeqn{\textrm{MDM} = \frac{1}{G} \sum_{g=1}^{G} d_{gM}}
#'
#'   Where \mjseqn{G} is the total number of entries in the cluster/group.
#'
#'   }
#'
#'   \subsection{Number of clusters}{\insertCite{diwan_core_1994}{SampleCore}
#'   proposed the number of clusters produced by a multivariate cluster
#'   analysis at a specific distance threshold as an estimate of the diversity.
#'
#'   }
#'
#' @template general-arg
#' @template log-arg
#' @template size-arg
#' @template dist-arg
#' @param method The allocation method. Either \code{"dist"} for constant or
#'   \code{"dist.prop"} for proportional or \code{"dist.log"} for logarithmic or
#'   \code{"dist.sqrt"} for square root allocation. See \strong{Methods}.
#' @param metric The metric to be computed from the distance matrix. Either
#'   \code{"mean"}, \code{"median"}, \code{"max"}, \code{"range"},
#'   \code{"mnnd"}, \code{"mdc"}, \code{"mdm"}, \code{"mstl"}, or
#'   \code{"nclust"}. See \strong{Metrics}.
#' @param clust.fun A function to generate clusters from a distance matrix and
#'   return the number of clusters.
#'
#' @template alloc-returns
#'
#' @importFrom vegan betadisper
#' @importFrom igraph E graph_from_adjacency_matrix mst
#' @importFrom stats as.dist setNames median
#' @export
#'
#' @references
#'
#' \insertAllCited
#'
#' @seealso \code{\link[SampleCore]{allocate.basic}},
#'   \code{\link[SampleCore]{allocate.diversity}}
#'
#' @examplesIf requireNamespace("cluster", quietly = TRUE) & requireNamespace("fastcluster", quietly = TRUE) & requireNamespace("dbscan", quietly = TRUE) & requireNamespace("biotools", quietly = TRUE)
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Prepare example data
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' library(cluster)
#'
#' # Get distance matrix
#' data("cassava_EC_gp")
#'
#' set.seed(123)
#' cassava_EC_gp <- cassava_EC_gp[sample(1:nrow(cassava_EC_gp), 500), ]
#'
#' quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW",
#'            "AVPW", "ARSR", "SRDM")
#' qual <- c("CUAL", "LNGS", "PTLC", "DSTA", "LFRT", "LBTEF", "CBTR", "NMLB",
#'           "ANGB", "CUAL9M", "LVC9M", "TNPR9M", "PL9M", "STRP", "STRC",
#'           "PSTR")
#'
#' data <- cassava_EC_gp
#'
#' # Convert qualitative data columns to factor
#' data[, qual] <- lapply(data[, qual], as.factor)
#'
#' # Standardise quantitative data column
#' data[, quant] <- lapply(data[, quant], function(x) {
#'   scale(x)[, 1]
#' })
#'
#' # Get the Gower's distance matrix
#' dist_matrix <- daisy(x = data[, c(qual, quant)],
#'                      metric = "gower")
#'
#' # Get data
#' data <- cassava_EC_gp
#' data <- cbind(genotypes = rownames(cassava_EC_gp), cassava_EC_gp)
#' row.names(data) <- NULL
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Custom clustering functions
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' # UPGMA with hclust
#' clust_fun_upgma <- function(x) {
#'   # Tree
#'   tree_out <- hclust(x, method = "average")
#'   # Clusters
#'   cutree(tree_out, h = 0.2)
#' }
#'
#' if (requireNamespace('fastcluster', quietly = TRUE)) {
#'   # Ward's minimum variance with fastcluster
#'   clust_fun_ward <- function(x) {
#'     # Tree
#'     tree_out <- fastcluster::hclust(x, method = "ward.D2")
#'     # Clusters
#'     cutree(tree_out, h = 0.2)
#'   }
#' }
#'
#' if (requireNamespace('dbscan', quietly = TRUE)) {
#'   # Density-based clustering with dbscan
#'   clust_fun_dbscan <- function(x) {
#'     clust_out <- dbscan::dbscan(x, eps = 0.25)
#'     # remove noise: TODO
#'     setNames(clust_out$cluster, labels(x))
#'   }
#' }
#'
#' if (requireNamespace('biotools', quietly = TRUE)) {
#'   # Tocher's sequential clustering
#'   clust_fun_tocher <- function(x) {
#'     clust_out <- biotools::tocher(x, algorithm = "sequential")
#'     setNames(clust_out$class, labels(x))
#'   }
#' }
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Diversity (Distance based) allocation
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' ## Mean
#' dist_out_mean <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "mean",
#'                     size = 0.2)
#' dist_out_mean
#'
#' ## Median
#' dist_out_median <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "median",
#'                     size = 0.2)
#' dist_out_median
#'
#' ## Maximum
#' dist_out_max <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "max",
#'                     size = 0.2)
#' dist_out_max
#'
#' ## Range
#' dist_out_range <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "range",
#'                     size = 0.2)
#' dist_out_range
#'
#' ## Mean nearest-neighbour distance
#' dist_out_mnnd <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "mnnd",
#'                     size = 0.2)
#' dist_out_mnnd
#'
#' ## Minimum spanning tree length
#' dist_out_mstl <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "mstl",
#'                     size = 0.2)
#' dist_out_mstl
#'
#' \donttest{
#'   ## Mean distance to centroid
#'   dist_out_mdc <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist",
#'                       dist.mat = dist_matrix, metric = "mdc",
#'                       size = 0.2)
#'   dist_out_mdc
#'
#'   ## Mean distance to median
#'   dist_out_mdm <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist",
#'                       dist.mat = dist_matrix, metric = "mdm",
#'                       size = 0.2)
#'   dist_out_mdm
#' }
#'
#' ## Number of clusters
#'
#' ### UPGMA with hclust
#' dist_out_nclust1 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_upgma,
#'                     size = 0.2)
#' dist_out_nclust1
#'
#' # Ward's minimum variance with fastcluster
#' dist_out_nclust2 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_ward,
#'                     size = 0.2)
#' dist_out_nclust2
#'
#'
#' # Density-based clustering with dbscan
#' dist_out_nclust3 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_dbscan,
#'                     size = 0.2)
#' dist_out_nclust3
#'
#'
#' # Tocher's sequential clustering
#' \donttest{
#'   dist_out_nclust4 <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist",
#'                       dist.mat = dist_matrix, metric = "nclust",
#'                       clust.fun = clust_fun_tocher,
#'                       size = 0.2)
#'   dist_out_nclust4
#' }
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Diversity (Distance based) & Proportional
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' ## Mean
#' dist_prop_out_mean <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "mean",
#'                     size = 0.2)
#' dist_prop_out_mean
#'
#' ## Median
#' dist_prop_out_median <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "median",
#'                     size = 0.2)
#' dist_prop_out_median
#'
#' ## Maximum
#' dist_prop_out_max <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "max",
#'                     size = 0.2)
#' dist_prop_out_max
#'
#' ## Range
#' dist_prop_out_range <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "range",
#'                     size = 0.2)
#' dist_prop_out_range
#'
#' ## Mean nearest-neighbour distance
#' dist_prop_out_mnnd <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "mnnd",
#'                     size = 0.2)
#' dist_prop_out_mnnd
#'
#' ## Minimum spanning tree length
#' dist_prop_out_mstl <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "mstl",
#'                     size = 0.2)
#' dist_prop_out_mstl
#'
#' \donttest{
#'   ## Mean distance to centroid
#'   dist_prop_out_mdc <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist.prop",
#'                       dist.mat = dist_matrix, metric = "mdc",
#'                       size = 0.2)
#'   dist_prop_out_mdc
#'
#'   ## Mean distance to median
#'   dist_prop_out_mdm <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist.prop",
#'                       dist.mat = dist_matrix, metric = "mdm",
#'                       size = 0.2)
#'   dist_prop_out_mdm
#' }
#'
#' ## Number of clusters
#'
#' ### UPGMA with hclust
#' dist_prop_out_nclust1 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_upgma,
#'                     size = 0.2)
#' dist_prop_out_nclust1
#'
#' # Ward's minimum variance with fastcluster
#' dist_prop_out_nclust2 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_ward,
#'                     size = 0.2)
#' dist_prop_out_nclust2
#'
#' # Density-based clustering with dbscan
#' dist_prop_out_nclust3 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.prop",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_dbscan,
#'                     size = 0.2)
#' dist_prop_out_nclust3
#'
#' # Tocher's sequential clustering
#' \donttest{
#'   dist_prop_out_nclust4 <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist.prop",
#'                       dist.mat = dist_matrix, metric = "nclust",
#'                       clust.fun = clust_fun_tocher,
#'                       size = 0.2)
#'   dist_prop_out_nclust4
#' }
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Diversity (Distance based) & Logarithmic
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' ## Mean
#' dist_log_out_mean <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "mean",
#'                     size = 0.2)
#' dist_log_out_mean
#'
#' ## Median
#' dist_log_out_median <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "median",
#'                     size = 0.2)
#' dist_log_out_median
#'
#' ## Maximum
#' dist_log_out_max <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "max",
#'                     size = 0.2)
#' dist_log_out_max
#'
#' ## Range
#' dist_log_out_range <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "range",
#'                     size = 0.2)
#' dist_log_out_range
#'
#' ## Mean nearest-neighbour distance
#' dist_log_out_mnnd <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "mnnd",
#'                     size = 0.2)
#' dist_log_out_mnnd
#'
#' ## Minimum spanning tree length
#' dist_log_out_mstl <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "mstl",
#'                     size = 0.2)
#' dist_log_out_mstl
#'
#' \donttest{
#'   ## Mean distance to centroid
#'   dist_log_out_mdc <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist.log",
#'                       dist.mat = dist_matrix, metric = "mdc",
#'                       size = 0.2)
#'   dist_log_out_mdc
#' }
#'
#' ## Mean distance to median
#' dist_log_out_mdm <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "mdm",
#'                     size = 0.2)
#' dist_log_out_mdm
#'
#' ## Number of clusters
#'
#' ### UPGMA with hclust
#' dist_log_out_nclust1 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_upgma,
#'                     size = 0.2)
#' dist_log_out_nclust1
#'
#' # Ward's minimum variance with fastcluster
#' dist_log_out_nclust2 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_ward,
#'                     size = 0.2)
#' dist_log_out_nclust2
#'
#' # Density-based clustering with dbscan
#' dist_log_out_nclust3 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.log",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_dbscan,
#'                     size = 0.2)
#' dist_log_out_nclust3
#'
#' # Tocher's sequential clustering
#' \donttest{
#'   dist_log_out_nclust4 <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist.log",
#'                       dist.mat = dist_matrix, metric = "nclust",
#'                       clust.fun = clust_fun_tocher,
#'                       size = 0.2)
#'   dist_log_out_nclust4
#' }
#'
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' # Diversity (Distance based) & Square root
#' #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' ## Mean
#' dist_sqrt_out_mean <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "mean",
#'                     size = 0.2)
#' dist_sqrt_out_mean
#'
#' ## Median
#' dist_sqrt_out_median <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "median",
#'                     size = 0.2)
#' dist_sqrt_out_median
#'
#' ## Maximum
#' dist_sqrt_out_max <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "max",
#'                     size = 0.2)
#' dist_sqrt_out_max
#'
#' ## Range
#' dist_sqrt_out_range <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "range",
#'                     size = 0.2)
#' dist_sqrt_out_range
#'
#' ## Mean nearest-neighbour distance
#' dist_sqrt_out_mnnd <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "mnnd",
#'                     size = 0.2)
#' dist_sqrt_out_mnnd
#'
#' ## Minimum spanning tree length
#' dist_sqrt_out_mstl <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "mstl",
#'                     size = 0.2)
#' dist_sqrt_out_mstl
#'
#' \donttest{
#'   ## Mean distance to centroid
#'   dist_sqrt_out_mdc <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist.sqrt",
#'                       dist.mat = dist_matrix, metric = "mdc",
#'                       size = 0.2)
#'   dist_sqrt_out_mdc
#'
#'   ## Mean distance to median
#'   dist_sqrt_out_mdm <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist.sqrt",
#'                       dist.mat = dist_matrix, metric = "mdm",
#'                       size = 0.2)
#'   dist_sqrt_out_mdm
#' }
#'
#' ## Number of clusters
#'
#' ### UPGMA with hclust
#' dist_sqrt_out_nclust1 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_upgma,
#'                     size = 0.2)
#' dist_sqrt_out_nclust1
#'
#' # Ward's minimum variance with fastcluster
#' dist_sqrt_out_nclust2 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_ward,
#'                     size = 0.2)
#' dist_sqrt_out_nclust2
#'
#' # Density-based clustering with dbscan
#' dist_sqrt_out_nclust3 <-
#'   allocate.distance(data = data, names = "genotypes",
#'                     group = "Cluster", method = "dist.sqrt",
#'                     dist.mat = dist_matrix, metric = "nclust",
#'                     clust.fun = clust_fun_dbscan,
#'                     size = 0.2)
#' dist_sqrt_out_nclust3
#'
#' # Tocher's sequential clustering
#' \donttest{
#'   dist_sqrt_out_nclust4 <-
#'     allocate.distance(data = data, names = "genotypes",
#'                       group = "Cluster", method = "dist.sqrt",
#'                       dist.mat = dist_matrix, metric = "nclust",
#'                       clust.fun = clust_fun_tocher,
#'                       size = 0.2)
#'   dist_sqrt_out_nclust4
#' }
#'
allocate.distance <- function(data, names, group,
                              dist.mat,
                              method = c("dist", "dist.prop",
                                         "dist.log", "dist.sqrt"),
                              metric = c("mean", "median", "max",
                                         "range", "mnnd", "mdc",
                                         "mdm", "mstl", "nclust"),
                              clust.fun = NULL,
                              log.base = exp(1),
                              size) {

  # Checks ----

  checks.sample.core(data = data, size = size,
                     names = names, group = group,
                     dist.mat = dist.mat,
                     qualitative = NULL,
                     log.base = log.base,
                     mode = "alloc")

  method <- match.arg(method)
  metric <- match.arg(metric)

  if (!is.null(clust.fun)) {
    if (!is.function(clust.fun)) {
      stop('"clust.fun" must be a function.')
    }

    clust_res <- clust.fun(dist.mat)
    if (!is.vector(clust_res) && !is.factor(clust_res)) {
      stop("clust.fun must return a vector or factor of cluster assignments.")
    }

    if (length(clust_res) != nrow(as.matrix(dist.mat))) {
      stop("clust.fun must return one cluster label per observation in x.")
    }
  }

  # Prepare data ----

  data[, group] <- droplevels(data[, group])

  dmat <- as.matrix(dist.mat)

  gp_memb <- data[, group]
  dist_labels <- labels(dist.mat)

  # mean_dist <- mean(dist.mat[upper.tri(dist.mat)], na.rm = TRUE)

  # Align labels
  names(gp_memb) <- dist_labels
  gp_memb <- gp_memb[rownames(dmat)]

  # Basic group stats ----

  gps <- levels(data[, group])
  gpsize <- summary(data[, group])
  gpcount <- length(gps)

  tcount <- nrow(data)

  # Compute group-wise distance metrics ----

  if (metric %in% c("mean", "median", "max",
                    "range", "mnnd", "mstl", "nclust") ) {

    group_dist_metric <- sapply(unique(gp_memb), function(g) {
      idx <- which(gp_memb == g)

      # extract submatrix for the group
      sub_d <- dmat[idx, idx]

      ## Mean/Average distance ----
      if (metric == "mean") {
        # take only upper triangle
        # avoid duplicates + diagonal
        out <- mean(sub_d[upper.tri(sub_d)], na.rm = TRUE)
      }

      ## Median distance ----
      if (metric == "median") {
        out <- median(sub_d[upper.tri(sub_d)])
      }

      ## Maximum distance ----
      if (metric == "max") {
        out <- max(sub_d[upper.tri(sub_d)])
      }

      ## Range ----
      if (metric == "range") {
        out <- diff(range(sub_d[upper.tri(sub_d)]))
      }

      ## Mean nearest-neighbour distance (MNND) ----
      # Captures clustering vs dispersion
      if (metric == "mnnd") {
        vals <- sub_d
        diag(vals) <- Inf # To exclude self-distances
        out <- mean(apply(vals, 1, min))
      }

      ## Minimum spanning tree (MSTL) length ----
      if (metric == "mstl") {
        # mst_mat <- ape::mst(sub_d)
        # mat <- as.matrix(sub_d)
        # sum(mst_mat * mat) / 2

        g <- igraph::graph_from_adjacency_matrix(
          as.matrix(sub_d),
          mode = "undirected",
          weighted = TRUE,
          diag = FALSE
        )

        mst_g <- igraph::mst(g)
        out <- sum(igraph::E(mst_g)$weight)
      }

      ## Diwan et al., 1994 (No. of clusters )----
      if (metric == "nclust") {
        # Handle groups with only 1 entry (cannot cluster)
        if (length(idx) == 1) {
          return(1)
        }

        # Return cluster membership
        clusters <- clust.fun(as.dist(sub_d))

        # count of unique clusters
        out <- length(unique(clusters))
      }

      return(out)

    })
  }

  ## Mean distance to centroid ----
  if (metric == "mdc") {
    # only for euclidean dist
    # coords <- cmdscale(sub_d, k = nrow(sub_d) - 1)
    # centroid <- colMeans(coords)
    # mean(sqrt(rowSums((coords - centroid) ^ 2)))

    bdout <- vegan::betadisper(d = dist.mat, group = data[, group],
                               type = "centroid")
    group_dist_metric <- bdout$group.distances
  }

  ## Mean distance to median ----
  if (metric == "mdm") {
    bdout <- vegan::betadisper(d = dist.mat, group = data[, group],
                               type = "median")
    group_dist_metric <- bdout$group.distances
  }

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
