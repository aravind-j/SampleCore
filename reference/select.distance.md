# Selection of Entries from Clusters/Groups on the basis of Genetic Distances

Select entries from cluster/groups in the entire collection by genetic
distance based sampling according to allocation specified.

## Usage

``` r
select.distance(
  data,
  names,
  group,
  alloc,
  dist.mat,
  always.selected,
  method = c("mean.medoid", "median.medoid", "nearest.centroid", "nearest.median",
    "mean.peripheral", "median.peripheral", "eccentricity", "farness.centrality",
    "kennard.stone", "duplex", "honigs", "farthest.point", "nearest.neighbour", "naes",
    "optim.medoid", "hclust.random", "hclust.medoid"),
  hclust.method = c("average", "single", "complete", "ward.D", "mcquitty", "median",
    "centroid", "ward.D2")
)
```

## Arguments

- data:

  The data as a data frame object. The data frame should possess one row
  per individual and columns with the individual names and multiple
  trait/character data.

- names:

  Name of column with the accession names as a character string.

- group:

  Name of column with the accession group/cluster names as a character
  string.

- alloc:

  A named numeric vector specifying the number of entries to be
  selected. Names should correspond to the levels of the "`"group"`
  column, and values indicate the number of elements to be selected from
  each level.

- dist.mat:

  A precomputed distance matrix of distance measures between the
  accessions in `data`.

- always.selected:

  Names of accessions to be always included in the core set as a
  character vector.

- method:

  The method for sampling accessions from each cluster/group. Either
  `"mean.medoid"`, `"median.medoid"`, `"nearest.centroid"`,
  `"nearest.median"`, `"mean.peripheral"`, `"median.peripheral"`,
  `"eccentricity"`, `"farness.centrality"`, `"kennard.stone"`,
  `"duplex"`, `"honigs"`, `"farthest.point"`, `"nearest.neighbour"`,
  `"naes"`, `"optim.medoid"`, `"hclust.random"` or `"hclust.medoid"`.
  See **Methods**.

- hclust.method:

  The hierarchical clustering method to be used. Either `"ward.D"`,
  `"ward.D2"`, `"single"`, `"complete"`, `"average"` (= UPGMA),
  `"mcquitty"` (= WPGMA), `"median"` (= WPGMC) or `"centroid"` (=
  UPGMC).

## Value

A named list where each element contains the selected entry identifiers
for a cluster/group.

## Details

For each cluster/group, entries are selected by several methods from
within-cluster/group genetic distances between accessions according to
the allocation provided (See **Methods**).

Entries listed as `always.selected` are mandatorily included in the
selection. Warnings are issued if requested allocation is smaller than
the number of always-selected entries in a cluster/group and/or when the
cluster/group does not contain enough remaining entries to fulfill the
allocation.

## Methods

### Centrality Based Methods

Selects accessions that are most representative/closest to the
cluster/group center.

#### Medoid-like Representative Sampling by Minimal Mean Distance

Selects medoid-like representatives as accessions with the smallest
average distance to all others within the group.

For each accession \\g\\, the mean distance to all other accessions
\\h\\ is computed as:

\\\bar{d}\_g = \frac{1}{G} \sum\_{h=1}^{G} d\_{gh}\\

Accessions are ranked by \\\bar{d}\_g\\ in ascending order and the top
\\n\\ are selected.

#### Medoid-like Representative Sampling by Minimal Median Distance

Selects medoid-like representatives as accessions with the smallest
median distance to all others within the group. This method is less
influenced by outliers.

For each accession \\g\\, the median distance to all other accessions
\\h\\ is computed as:

\\\tilde{d}\_g = \text{median}\_{h=1,\dots,G}(d\_{gh})\\

Accessions are ranked by \\\tilde{d}\_g\\ in ascending order and the top
\\n\\ are selected.

#### Representative Sampling by Proximity to Group Centroid

Selects accessions closest to the group centroid in principal coordinate
space, computed via multivariate dispersion analysis using
[`betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html).

The distance of each accession \\g\\ to the group centroid \\C\\ in PCoA
space is:

\\\delta_g = \\ \mathbf{p}\_g - \mathbf{c} \\\\

Where \\\mathbf{p}\_g\\ is the PCoA coordinate vector of accession \\g\\
and \\\mathbf{c}\\ is the group centroid. Accessions are ranked by
\\\delta_g\\ in ascending order and the top \\n\\ are selected.

#### Representative Sampling by Proximity to Group Spatial Median

Selects accessions closest to the group spatial median in principal
coordinate space, computed via multivariate dispersion analysis using
[`betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html).

The distance of each accession \\g\\ to the group spatial median \\M\\
is:

\\\delta_g^\* = \\ \mathbf{p}\_g - \mathbf{m} \\\\

where \\\mathbf{m}\\ is the spatial median of the group in PCoA space.
Accessions are ranked by \\\delta_g^\*\\ in ascending order and the top
\\n\\ are selected.

### Peripheral/Extremity Based Methods

Selects accessions that are most dissimilar from the rest in a
cluster/group i.e. the accessions which are in the boundary or outliers.

#### Peripheral Sampling by Maximal Mean Distance

Selects the most peripheral accessions as those with the largest average
distance to all others within the group.

\\\bar{d}\_g = \frac{1}{G} \sum\_{h=1}^{G} d\_{gh}\\

Accessions are ranked by \\\bar{d}\_g\\ in descending order and the top
\\n\\ are selected.

#### Peripheral Sampling by Maximal Median Distance

Selects the most peripheral accessions as those with the largest median
distance to all others within the group.

\\\tilde{d}\_g = \text{median}\_{h=1,\dots,G}(d\_{gh})\\

Accessions are ranked by \\\tilde{d}\_g\\ in descending order and the
top \\n\\ are selected.

#### Peripheral Sampling by Maximal Eccentricity

Selects accessions with the largest eccentricity — the maximum distance
to any other accession in the group.

\\e_g = \max\_{h=1,\dots,G} d\_{gh}\\

Accessions are ranked by \\e_g\\ in descending order and the top \\n\\
are selected. Eccentricity captures the worst-case dissimilarity of an
accession rather than its average behaviour.

#### Peripheral Sampling by Maximal Farness Centrality

Selects accessions with the greatest total distance to all others, i.e.
those most remote from the rest of the group.

\\f_g = \sum\_{h=1}^{G} d\_{gh}\\

Accessions are ranked by \\f_g\\ in descending order and the top \\n\\
are selected. Farness centrality is proportional to \\\bar{d}\_g\\ and
differs from `mean.peripheral` only in that it uses the raw sum rather
than the mean, producing identical rankings.

### Space-Filling/Coverage Methods

Select accessions that are spread maximally across the feature space in
a cluster/group i.e. diversity sampling.

#### Space-Filling Sampling via the Kennard-Stone Algorithm

Selects \\n\\ accessions that maximally and uniformly cover the distance
space via the Kennard-Stone algorithm (See
[`kenStone`](https://rdrr.io/pkg/prospectr/man/kenStone.html)).

Starting from the pair of accessions with the largest pairwise distance:

\\\lbrace g_1, g_2 \rbrace = \underset{g,h}{\arg\max}\\ d\_{gh}\\

each subsequent accession \\g_k\\ is selected by maximising its minimum
distance to the already-selected set \\S\\:

\\g_k = \underset{g \notin S}{\arg\max} \min\_{s \in S} d\_{gs}\\

This greedy procedure ensures even space coverage without relying on
cluster structure.

#### Space-Filling Sampling via the DUPLEX Algorithm

Extends the Kennard-Stone algorithm to simultaneously construct a model
set and a test set with similar distributions
([duplex](https://rdrr.io/pkg/prospectr/man/duplex.html)). Accessions
are selected using Mahalanobis distance:

\\d_M(g, h) = \sqrt{(\mathbf{x}\_g - \mathbf{x}\_h)^\top \Sigma^{-1}
(\mathbf{x}\_g - \mathbf{x}\_h)}\\

where \\\Sigma\\ is the covariance matrix. At each step, the pair
maximising \\d_M\\ is split alternately between model and test sets,
ensuring both sets span the full feature space.

#### Space-Filling Sampling via the Honigs Algorithm

Selects \\n\\ accessions sequentially by maximising dissimilarity to the
already-selected set
([honigs](https://rdrr.io/pkg/prospectr/man/honigs.html))

At each step \\k\\, the accession \\g_k\\ maximising total distance to
all previously selected accessions \\S\\ is chosen:

\\g_k = \underset{g \notin S}{\arg\max} \sum\_{s \in S} d\_{gs}\\

This favours accessions that are collectively most dissimilar to the
current selection, producing broad coverage of the distance space.

.

#### Space-Filling Sampling via Farthest-Point (Max-Min) Algorithm

Selects \\n\\ accessions by iteratively maximising the minimum distance
to the current selected set — also known as the max-min or
farthest-point sampling algorithm.

\\g_k = \underset{g \notin S}{\arg\max} \min\_{s \in S} d\_{gs}\\

This is equivalent to Kennard-Stone but without the symmetric
initialisation step. It provides a deterministic, greedy approximation
to the \\k\\-centre problem:

\\\min\_{S \subset G,\\ \|S\|=n} \max\_{g \in G} \min\_{s \in S}
d\_{gs}\\

### Density Based Methods

Select points based on local neighbourhood density.

#### Density-Based Sampling by Minimal Nearest-Neighbour Distance

Selects accessions residing in the densest regions of the distance
space, identified as those with the smallest nearest-neighbour distance.

For each accession \\g\\, the nearest-neighbour distance is:

\\\text{nn}\_g = \min\_{h \neq g} d\_{gh}\\

Accessions are ranked by \\\text{nn}\_g\\ in ascending order and the top
\\n\\ are selected. Small \\\text{nn}\_g\\ indicates that \\g\\ resides
in a dense cluster; this method preferentially samples from high-density
regions.

### Cluster Based Methods

These methods partition the cluster/group space into
sub-clusters/groups, then samples from each one.

#### Globally Optimal Medoid Sampling via Partitioning Around Medoids (PAM)

Selects a set of \\n\\ medoids that jointly minimise the total distance
of every accession to its nearest medoid, via
[`pam`](https://rdrr.io/pkg/cluster/man/pam.html).

The objective function minimised is:

\\\min\_{S \subset G,\\ \|S\|=n} \sum\_{g=1}^{G} \min\_{s \in S}
d\_{gs}\\

Unlike `"mean.medoid"`, medoids are co-optimised as a set, ensuring they
collectively represent the full distribution of the group rather than
independently scoring each accession.

#### Cluster-Based Sampling via K-means (Naes Method)

Partitions accessions into \\n\\ clusters via k-means applied to the
distance matrix (See
[`naes`](https://rdrr.io/pkg/prospectr/man/naes.html)), then selects the
accession closest to each cluster centre as the representative.

The k-means objective minimised is:

\\\min \sum\_{k=1}^{n} \sum\_{g \in C_k} d\_{g, \mu_k}^2\\

where \\C_k\\ is the \\k\\-th cluster and \\\mu_k\\ is its centre. One
representative per cluster is returned, ensuring broad, partition-aware
coverage.

#### Cluster-Based Sampling via Hierarchical Clustering with Random Selection

Partitions accessions into \\n\\ clusters by cutting a hierarchical
clustering dendrogram at height \\k = n\\, then randomly samples one
accession from each cluster.

The dendrogram is built by agglomerative hierarchical clustering using
the linkage criterion specified by
[`hclust`](https://rdrr.io/r/stats/hclust.html). For clusters \\C_1,
\dots, C_n\\, one accession is drawn uniformly at random from each:

\\g_k \sim \text{Uniform}(C_k), \quad k = 1, \dots, n\\

This introduces stochasticity within a structured partition, balancing
coverage with randomness.

#### Cluster-Based Sampling via Hierarchical Clustering with Medoid Selection

Partitions accessions into \\n\\ clusters by cutting a hierarchical
clustering dendrogram at height \\k = n\\, then selects the
within-cluster medoid as the representative of each cluster.

For each cluster \\C_k\\, the medoid is the accession minimising total
within-cluster distance:

\\g_k^\* = \underset{g \in C_k}{\arg\min} \sum\_{h \in C_k} d\_{gh}\\

This combines the structured partitioning of hierarchical clustering
with deterministic, centrality-based representative selection.
