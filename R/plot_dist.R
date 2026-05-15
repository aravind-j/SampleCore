#' Plot a distance matrix as a 2D projection
#'
#' Reduces a distance matrix to two dimensions using Classical MDS, Isotonic
#' MDS, or t-SNE, and returns a \code{ggplot2} scatter plot in which proximity
#' reflects similarity. Points can optionally be highlighted or split into facet
#' panels by group.
#'
#' @param d A distance matrix of class \code{\link[stats]{dist}}. Labels must be
#'   set (i.e. \code{labels(d)} must not be \code{NULL}). Duplicate labels are
#'   not permitted.
#' @param method Character string specifying the dimensionality-reduction
#'   method. One of:
#'   \describe{
#'     \item{\code{"cmds"}}{Classical (metric) Multidimensional Scaling via
#'       \code{\link[stats]{cmdscale}}. This is the default.}
#'     \item{\code{"isomds"}}{Non-metric (isotonic) MDS via
#'       \code{\link[MASS]{isoMDS}}. Automatically falls back to
#'       \code{"cmds"} with a message when \code{n < 3}.}
#'     \item{\code{"tsne"}}{t-distributed Stochastic Neighbour Embedding via
#'       \code{\link[Rtsne]{Rtsne}}. Perplexity is set automatically to
#'       \code{min(30, floor((n - 1) / 3))}.}
#'   }
#' @param highlight Optional character vector of labels to highlight in the
#'   plot. Matching identifiers are plotted in \strong{red}; all others in
#'   black. \code{NULL} (default) disables highlighting. Every value must be
#'   present in \code{labels(d)}.
#' @param gp Optional named character vector mapping labels to group names
#'   (\code{names(gp)} = labels, values = group names). When supplied, the plot
#'   is split into one facet panel per group via
#'   \code{\link[ggplot2]{facet_wrap}}. The set of names must match
#'   \code{labels(d)} exactly. \code{NULL} (default) produces a single panel.
#' @param point.alpha Alpha transparency value for points.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object. The plot can be further
#'   customised with standard \pkg{ggplot2} additions before printing or saving.
#'
#' @seealso \code{\link[stats]{cmdscale}}, \code{\link[MASS]{isoMDS}},
#' \code{\link[Rtsne]{Rtsne}}, \code{\link[ggplot2]{ggplot}}
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_bw theme element_rect facet_wrap
#' @importFrom MASS isoMDS
#' @importFrom Rtsne Rtsne
#' @export
#'
#' @examples
#' # Basic usage with the built-in eurodist dataset
#' plot_dist(eurodist)
#'
#' # Non-metric MDS with two highlighted cities
#' plot_dist(eurodist, method = "isomds",
#'           highlight = c("Madrid", "Rome"))
#'
#' # Classical MDS split by a user-defined grouping
#' regions <-
#'   c(Athens = "South",  Barcelona = "South", Brussels = "North",
#'     Calais = "North",  Cherbourg = "North", Cologne = "North",
#'     Copenhagen = "North", Geneva = "South", Gibraltar = "South",
#'     Hamburg = "North", `Hook of Holland` = "North", Lisbon = "South",
#'     Lyons = "South",  Madrid = "South",  Marseilles = "South",
#'     Milan = "South",  Munich = "North",  Paris = "North",
#'     Rome = "South",   Stockholm = "North", Vienna = "North")
#'
#' plot_dist(eurodist, method = "cmds", gp = regions,
#'           highlight = c("Madrid", "Cherbourg", "Rome", "Brussels"))
#'
plot_dist <- function(d,
                      method = c("cmds", "isomds", "tsne"),
                      highlight = NULL,
                      gp = NULL, point.alpha = 0.8) {

  # Checks ----
  method <- match.arg(method)

  d_labs <- attr(d, "Labels")

  if (!is.null(gp)) {
    missing <- setdiff(d_labs, names(gp))
    if (length(missing) > 0) {
      warning('Some objects in distance matrix "d" have no group assignment: ',
              paste(missing, collapse = ", "))
    }
  }

  # check if 'd' is a distance matrix
  if (!inherits(d, "dist")) {
    stop('"d" should be a distance matrix of class "dist".')
  }
  dist_labels <- d_labs # No dups as it is a distance matrix

  if (is.null(dist_labels)) {
    stop('Labels are missing in distance matrix "d".')
  }

  if (!is.null(highlight)) {
    # check if 'highlight' is a character vector
    if (!is.character(highlight)) {
      stop('"highlight" should be a character vector.')
    }

    # check if highlight is present in the entire set
    if (any(!(highlight %in% d_labs))) {
      alsel_miss <- highlight[!(highlight %in% d_labs)]
      stop(paste('The following entry/entries specified in "highlight" ',
                 'are not present in "data":\n',
                 paste(alsel_miss, collapse = ", "),
                 sep = ""))
    }
  }

  # check point.alpha argument is numeric vector of unit length
  if (!(is.numeric(point.alpha) && length(point.alpha) == 1)) {
    stop('"point.alpha" should be a numeric vector of unit length.')
  }

  n <- attr(d, "Size")

  # fallback for small matrices
  if (method == "isomds" && n < 3) {
    message("isoMDS requires n >= 3; falling back to cmds.")
    method <- "cmds"
  }

  if (method == "tsne" && n < 4) {
    message("tsne requires n >= 4; falling back to cmds.")
    method <- "cmds"
  }

  # Dimensional reduction ----
  if (method == "cmds") {
    fit <- cmdscale(d, k = 2)
  }

  if (method == "isomds") {
    fit <- MASS::isoMDS(d, trace = FALSE)$points
  }

  if (method == "tsne") {

    prplx <- max(1, min(30, floor((n - 1) / 3)))

    fit <- Rtsne::Rtsne(as.matrix(d),
      is_distance = TRUE, perplexity = prplx)$Y

    rownames(fit) <- d_labs
  }

  # Labels ----
  labs <- d_labs

  # Plotting dataframe ----
  df <- data.frame(sample = labs,
                   X = fit[, 1],
                   Y = fit[, 2])

  # Highlighting ----
  df$highlight <- if (!is.null(highlight) &&
                      length(highlight) > 0) {
    ifelse(df$sample %in% highlight, "highlight", "other")
  } else {
    "other"
  }

  # Ensure highlighted points are plotted last
  df$highlight <- factor(df$highlight,
                         levels = c("other", "highlight"))

  df <- df[order(df$highlight), ]

  # Groups ----
  if (!is.null(gp)) {

    # gp should be a named vector:
    # names(gp) = sample labels
    # values(gp) = group labels

    if (is.null(names(gp))) {
      stop("gp must be a named vector with sample labels as names.")
    }

    df$group <- gp[df$sample]

  } else {
    df$group <- "all"
  }

  # ggplot ----

  axis_label <- switch(method, cmds = "PCoA", isomds = "MDS", tsne = "t-SNE")

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = X, y = Y,
                                    color = highlight)) +
    ggplot2::geom_point(size = 3, alpha = point.alpha) +
    ggplot2::scale_color_manual(values = c(other = "black",
        highlight = "red")) +
    ggplot2::labs(
      x = paste0(axis_label, " 1"), y = paste0(axis_label, " 2"),
      color = NULL,
      title = paste("Distance projection:", method)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "grey95"))

  # Faceting ----
  if (!is.null(gp)) {
    p <- p +
      ggplot2::facet_wrap(~group)
  }

  return(p)
}

