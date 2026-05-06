#' Allocation of Accessions on the Basis of Within Cluster/Group Distance-based
#' Diversity
#'
#' Estimate the number of accessions to be allocated from each cluster/group in
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
#'   the basis of diversity of accessions within that cluster/group alone or in
#'   combination with the size of the cluster/group (See
#'   \strong{\code{Methods}}).
#'
#'   The within-cluster/group diversity is estimated as several metrics from the
#'   within cluster/group genetic distances between accessions (See
#'   \strong{\code{Metrics}}).
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
#'   maximum and range of genetic distances between accessions in a cluster.}
#'
#'   \subsection{Mean nearest-neighbour distance (\mjseqn{MNND})}{It is the
#'   average, across all accessions, of the distance to each accession’s
#'   closest other accession (\mjseqn{d_{g_{min}}}), based on a genetic given
#'   distance matrix \insertCite{clark_distance_1954}{SampleCore}.
#'
#'   For each accession, the nearest-neighbour distance (\mjseqn{d_{g_{min}}})
#'   is the smallest non-zero distance with any other accession.
#'
#'   \mjsdeqn{d_{g_{min}} = \min_{h \ne g} d_{gh}}
#'
#'   The Mean nearest-neighbour distance (\mjseqn{MNND}) can then be computed
#'   as:
#'
#'   \mjsdeqn{\textrm{MNND} = \frac{1}{G} \sum_{g=1}^{G} d_g}
#'
#'   Where, (\mjseqn{g}) is the index of an accession in a genetic distance
#'   matrix, \mjseqn{h} is the index of all other genotypes and \mjseqn{G} is
#'   the total number of genotypes in a cluster/group.
#'
#'   }
#'
#'   \subsection{Minimum spanning tree length (\mjseqn{MSTL})}{It is defined as
#'   the sum of edge weights in the minimum spanning tree constructed from the
#'   genetic distance matrix of accessions within a cluster/group.  A minimum
#'   spanning tree (MST) connects all accessions such that the total distance
#'   is minimized and no cycles are formed. It represents the most efficient
#'   way to connect all accessions based on pairwise genetic distances
#'   \insertCite{gower_minimum_1969}{SampleCore}.
#'
#'   For genetic distance \mjseqn{d_{gh}} between accessions \mjseqn{g} and
#'   \mjseqn{h}, the MST is a subset of edges that connects all
#'   \mjseqn{G} accessions with exactly \mjseqn{G - 1} edges and minimum total
#'   weight. The MST length (\mjseqn{MSTL}) can then be computed as:
#'
#'   \mjsdeqn{\textrm{MSTL} = \sum_{(g,h) \in \mathcal{T}} d_{gh}}
#'
#'   Where \mjseqn{\mathcal{T}} denotes the set of edges in the MST.
#'
#'   }
#'
#'   \subsection{Mean distance to centroid and median (\mjseqn{MDC},
#'   \mjseqn{MDM})}{These quantify the average dispersion of accessions within a
#'   cluster/group relative to a central point in multivariate space derived
#'   from the genetic distance matrix.
#'
#'   The centroid represents the multivariate mean position of all accessions
#'   in a cluster
#'   \insertCite{sokal_principles_1963,sneath_numerical_1973}{SampleCore}.,
#'   whereas the median (spatial median) provides a robust central location
#'   that is less influenced by extreme values
#'   \insertCite{bradley_constrained_1999}{SampleCore}.
#'
#'   For \mjseqn{d_{gC}} and \mjseqn{d_{gM}} distances of accession
#'   \mjseqn{g} from the centroid \mjseqn{C} and median \mjseqn{M},
#'   respectively. These measures are computed as:
#'
#'   \mjsdeqn{\textrm{MDC} = \frac{1}{G} \sum_{g=1}^{G} d_{gC}}
#'
#'   \mjsdeqn{\textrm{MDM} = \frac{1}{G} \sum_{g=1}^{G} d_{gM}}
#'
#'   Where \mjseqn{G} is the total number of accessions in the cluster/group.
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
#' @param method The allocation method. Either \code{"dist"} for constant or
#'   \code{"dist.prop"} for proportional or \code{"dist.log"} for logarithmic or
#'   \code{"dist.sqrt"} for square root allocation.
#' @param clust.fun A function to generate clusters from a distance matrix and
#'   return the number of clusters.
#'
#' @returns
#'
#' @importFrom vegan betadisper
#' @importFrom igraph E graph_from_adjacency_matrix mst
#' @export
#'
#' @references
#'
#' \insertAllCited
#'
#' @examples
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
                     quantitative = NULL,
                     qualitative = NULL,
                     log.base = log.base)

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
                    "range", "mnnd", "mstl") ) {

    group_dist_metric <- sapply(unique(gp_memb), function(g) {
      idx <- which(gp_memb == g)

      # extract submatrix for the group
      sub_d <- dmat[idx, idx]

      ## Mean/Average distance ----
      if (metric == "mean") {
        # take only upper triangle
        # avoid duplicates + diagonal
        mean(sub_d[upper.tri(sub_d)], na.rm = TRUE)
      }

      ## Median distance ----
      if (metric == "median") {
        median(sub_d[upper.tri(sub_d)])
      }

      ## Maximum distance ----
      if (metric == "max") {
        max(sub_d[upper.tri(sub_d)])
      }

      ## Range ----
      if (metric == "range") {
        diff(range(sub_d[upper.tri(sub_d)]))
      }

      ## Mean nearest-neighbour distance (MNND) ----
      # Captures clustering vs dispersion
      if (metric == "mnnd") {
        vals <- sub_d
        diag(vals) <- Inf # To exclude self-distances
        mean(apply(vals, 1, min))
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
        sum(igraph::E(mst_g)$weight)
      }

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

  ## Diwan et al., 1994 (No. of clusters )----
  if (metric == nclust) {
    # Handle groups with only 1 accession (cannot cluster)
    if (length(idx) == 1) {
      return(1)
    }

    # Return cluster membership
    clusters <- clust.fun(as.dist(sub_d))

    # count of unique clusters
    length(unique(clusters))
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
