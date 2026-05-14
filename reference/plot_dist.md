# Plot a distance matrix as a 2D projection

Reduces a distance matrix to two dimensions using Classical MDS,
Isotonic MDS, or t-SNE, and returns a `ggplot2` scatter plot in which
proximity reflects similarity. Points can optionally be highlighted or
split into facet panels by group.

## Usage

``` r
plot_dist(
  d,
  method = c("cmds", "isomds", "tsne"),
  highlight = NULL,
  gp = NULL,
  point.alpha = 0.6
)
```

## Arguments

- d:

  A distance matrix of class
  [`dist`](https://rdrr.io/r/stats/dist.html). Labels must be set (i.e.
  `labels(d)` must not be `NULL`). Duplicate labels are not permitted.

- method:

  Character string specifying the dimensionality-reduction method. One
  of:

  `"cmds"`

  :   Classical (metric) Multidimensional Scaling via
      [`cmdscale`](https://rdrr.io/r/stats/cmdscale.html). This is the
      default.

  `"isomds"`

  :   Non-metric (isotonic) MDS via
      [`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html).
      Automatically falls back to `"cmds"` with a message when `n < 3`.

  `"tsne"`

  :   t-distributed Stochastic Neighbour Embedding via
      [`Rtsne`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html). Perplexity is
      set automatically to `min(30, floor((n - 1) / 3))`.

- highlight:

  Optional character vector of labels to highlight in the plot. Matching
  identifiers are plotted in **red**; all others in black. `NULL`
  (default) disables highlighting. Every value must be present in
  `labels(d)`.

- gp:

  Optional named character vector mapping labels to group names
  (`names(gp)` = labels, values = group names). When supplied, the plot
  is split into one facet panel per group via
  [`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).
  The set of names must match `labels(d)` exactly. `NULL` (default)
  produces a single panel.

- point.alpha:

  Alpha transparency value for points.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object. The plot can be further customised with standard ggplot2
additions before printing or saving.

## See also

[`cmdscale`](https://rdrr.io/r/stats/cmdscale.html),
[`isoMDS`](https://rdrr.io/pkg/MASS/man/isoMDS.html),
[`Rtsne`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html),
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Examples

``` r
# Basic usage with the built-in eurodist dataset
plot_dist(eurodist)


# Non-metric MDS with two highlighted cities
plot_dist(eurodist, method = "isomds",
          highlight = c("Madrid", "Rome"))


# Classical MDS split by a user-defined grouping
regions <-
  c(Athens = "South",  Barcelona = "South", Brussels = "North",
    Calais = "North",  Cherbourg = "North", Cologne = "North",
    Copenhagen = "North", Geneva = "South", Gibraltar = "South",
    Hamburg = "North", `Hook of Holland` = "North", Lisbon = "South",
    Lyons = "South",  Madrid = "South",  Marseilles = "South",
    Milan = "South",  Munich = "North",  Paris = "North",
    Rome = "South",   Stockholm = "North", Vienna = "North")

plot_dist(eurodist, method = "cmds", gp = regions,
          highlight = c("Madrid", "North", "Rome", "Brussels"))
#> Error in plot_dist(eurodist, method = "cmds", gp = regions, highlight = c("Madrid",     "North", "Rome", "Brussels")): The following entry/entries specified in "highlight" are not present in "data":
#> North
```
