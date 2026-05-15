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
  always.selected = NULL,
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
average distance to all others within the group (Kaufman and Rousseeuw
1987; Kaufman and Rousseeuw 1990) .

For each accession \\g\\, the mean distance to all other accessions
\\h\\ is computed as:

\\\bar{d}\_g = \frac{1}{G} \sum\_{h=1}^{G} d\_{gh}\\

Accessions are ranked by \\\bar{d}\_g\\ in ascending order and the top
\\n\\ are selected.

#### Medoid-like Representative Sampling by Minimal Median Distance

Selects medoid-like representatives as accessions with the smallest
median distance to all others within the group. This method is less
influenced by outliers (Kaufman and Rousseeuw 1987; Kaufman and
Rousseeuw 1990) .

For each accession \\g\\, the median distance to all other accessions
\\h\\ is computed as:

\\\tilde{d}\_g = \text{median}\_{h=1,\dots,G}(d\_{gh})\\

Accessions are ranked by \\\tilde{d}\_g\\ in ascending order and the top
\\n\\ are selected.

#### Representative Sampling by Proximity to Group Centroid

Selects accessions closest to the group centroid in principal coordinate
space, computed via multivariate dispersion analysis using
[`betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html)
(Anderson 2006; Anderson et al. 2006) .

The distance of each accession \\g\\ to the group centroid \\C\\ in PCoA
space is:

\\\delta_g = \\ \mathbf{p}\_g - \mathbf{c} \\\\

Where \\\mathbf{p}\_g\\ is the PCoA coordinate vector of accession \\g\\
and \\\mathbf{c}\\ is the group centroid. Accessions are ranked by
\\\delta_g\\ in ascending order and the top \\n\\ are selected.

#### Representative Sampling by Proximity to Group Spatial Median

Selects accessions closest to the group spatial median in principal
coordinate space, computed via multivariate dispersion analysis using
[`betadisper`](https://vegandevs.github.io/vegan/reference/betadisper.html)
(O'Neill and Mathews 2000) .

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
distance to all others within the group (Kaufman and Rousseeuw 1987;
Kaufman and Rousseeuw 1990) .

\\\bar{d}\_g = \frac{1}{G} \sum\_{h=1}^{G} d\_{gh}\\

Accessions are ranked by \\\bar{d}\_g\\ in descending order and the top
\\n\\ are selected.

#### Peripheral Sampling by Maximal Median Distance

Selects the most peripheral accessions as those with the largest median
distance to all others within the group (Kaufman and Rousseeuw 1987;
Kaufman and Rousseeuw 1990) .

\\\tilde{d}\_g = \text{median}\_{h=1,\dots,G}(d\_{gh})\\

Accessions are ranked by \\\tilde{d}\_g\\ in descending order and the
top \\n\\ are selected.

#### Peripheral Sampling by Maximal Eccentricity

Selects accessions with the largest eccentricity — the maximum distance
to any other accession in the group (Hage and Harary 1995) .

\\e_g = \max\_{h=1,\dots,G} d\_{gh}\\

Accessions are ranked by \\e_g\\ in descending order and the top \\n\\
are selected. Eccentricity captures the worst-case dissimilarity of an
accession rather than its average behaviour.

#### Peripheral Sampling by Maximal Farness Centrality

Selects accessions with the greatest total distance to all others, i.e.
those most remote from the rest of the group (Sabidussi 1966) .

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
space via the Kennard-Stone algorithm (Kennard and Stone 1969) (See
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
set and a test set with similar distributions (Kennard and Stone 1969;
Snee 1977) ([duplex](https://rdrr.io/pkg/prospectr/man/duplex.html)).
Accessions are selected using Mahalanobis distance:

\\d_M(g, h) = \sqrt{(\mathbf{x}\_g - \mathbf{x}\_h)^\top \Sigma^{-1}
(\mathbf{x}\_g - \mathbf{x}\_h)}\\

where \\\Sigma\\ is the covariance matrix. At each step, the pair
maximising \\d_M\\ is split alternately between model and test sets,
ensuring both sets span the full feature space.

#### Space-Filling Sampling via the Honigs Algorithm

Selects \\n\\ accessions sequentially by maximising dissimilarity to the
already-selected set (Honigs et al. 1985)
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
farthest-point sampling algorithm (Gonzalez 1985; Dyer and Frieze 1985;
Hochbaum and Shmoys 1985) .

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
space, identified as those with the smallest nearest-neighbour distance
(Cover and Hart 1967; Fix and Hodges 1989) .

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
accession closest to each cluster centre as the representative (Naes
1987; Naes et al. 2017) .

The k-means objective minimised is:

\\\min \sum\_{k=1}^{n} \sum\_{g \in C_k} d\_{g, \mu_k}^2\\

where \\C_k\\ is the \\k\\-th cluster and \\\mu_k\\ is its centre. One
representative per cluster is returned, ensuring broad, partition-aware
coverage.

#### Cluster-Based Sampling via Hierarchical Clustering with Random Selection

Partitions accessions into \\n\\ clusters by cutting a hierarchical
clustering dendrogram at height \\k = n\\, then randomly samples one
accession from each cluster (Ward 1963; Li et al. 2002) .

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
within-cluster medoid as the representative of each cluster (Kaufman and
Rousseeuw 1987; Ward 1963) .

For each cluster \\C_k\\, the medoid is the accession minimising total
within-cluster distance:

\\g_k^\* = \underset{g \in C_k}{\arg\min} \sum\_{h \in C_k} d\_{gh}\\

This combines the structured partitioning of hierarchical clustering
with deterministic, centrality-based representative selection.

## References

Anderson MJ (2006). “Distance-based tests for homogeneity of
multivariate dispersions.” *Biometrics*, **62**(1), 245–253.  
  
Anderson MJ, Ellingsen KE, McArdle BH (2006). “Multivariate dispersion
as a measure of beta diversity.” *Ecology Letters*, **9**(6), 683–693.  
  
Cover T, Hart P (1967). “Nearest neighbor pattern classification.” *IEEE
Transactions on Information Theory*, **13**(1), 21–27.  
  
Dyer ME, Frieze AM (1985). “A simple heuristic for the *p*-centre
problem.” *Operations Research Letters*, **3**(6), 285–288.  
  
Fix E, Hodges JL (1989). “Discriminatory analysis - Nonparametric
discrimination: Consistency properties.” *International Statistical
Review / Revue Internationale de Statistique*, **57**(3), 238–247.  
  
Gonzalez TF (1985). “Clustering to minimize the maximum intercluster
distance.” *Theoretical Computer Science*, **38**, 293–306.  
  
Hage P, Harary F (1995). “Eccentricity and centrality in networks.”
*Social Networks*, **17**(1), 57–63.  
  
Hochbaum DS, Shmoys DB (1985). “A best possible heuristic for the
*K*-center problem.” *Mathematics of Operations Research*, **10**(2),
180–184.  
  
Honigs DE, Hieftje GM, Mark HL, Hirschfeld TB (1985). “Unique-sample
selection via near-infrared spectral subtraction.” *Analytical
Chemistry*, **57**(12), 2299–2303.  
  
Kaufman L, Rousseeuw PJ (1990). *Finding Groups in Data: An Introduction
to Cluster Analysis*, Wiley Series in Probability and Statistics, 1
edition. Wiley. ISBN 978-0-471-87876-6 978-0-470-31680-1.  
  
Kaufman P, Rousseeuw PJ (1987). “Clustering by means of medoids.” In
Dodge Y (ed.), *Proceedings of the Statistical Data Analysis Based on
the L1 Norm Conference, Neuchatel, Switzerland*, volume 31, 405–416.  
  
Kennard RW, Stone LA (1969). “Computer aided design of experiments.”
*Technometrics*, **11**(1), 137–148.  
  
Li Z, Zhang H, Zeng Y, Yang Z, Shen S, Sun C, Wang X (2002). “Studies on
sampling schemes for the establishment of core collection of rice
landraces in Yunnan, China.” *Genetic Resources and Crop Evolution*,
**49**(1), 67–74.  
  
Naes T (1987). “The design of calibration in near infra-red reflectance
analysis by clustering.” *Journal of Chemometrics*, **1**(2), 121–134.  
  
Naes T, Isaksson T, Fearn T, Davies T (2017). *A User-Friendly Guide to
Multivariate Calibration and Classification*, Second edition edition. IM
Publications LLP, Chichester. ISBN 978-1-906715-25-0.  
  
O'Neill ME, Mathews K (2000). “A weighted least squares approach to
levene's test of homogeneity of variance.” *Australian & New Zealand
Journal of Statistics*, **42**(1), 81–100.  
  
Sabidussi G (1966). “The centrality index of a graph.” *Psychometrika*,
**31**(4), 581–603.  
  
Snee RD (1977). “Validation of regression models: Methods and examples.”
*Technometrics*, **19**(4), 415–428.  
  
Ward JH (1963). “Hierarchical grouping to optimize an objective
function.” *Journal of the American Statistical Association*,
**58**(301), 236–244.

## See also

[`select.random`](https://aravind-j.github.io/SampleCore/reference/select.random.md),
[`select.diversity`](https://aravind-j.github.io/SampleCore/reference/select.diversity.md)

## Examples

``` r

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare example data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(cluster)
library(ggplot2)

data(cassava_EC_gp)

data <- cassava_EC_gp

quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW", "AVPW",
           "ARSR", "SRDM")
qual <- c("CUAL", "LNGS", "PTLC", "DSTA", "LFRT", "LBTEF", "CBTR", "NMLB",
          "ANGB", "CUAL9M", "LVC9M", "TNPR9M", "PL9M", "STRP", "STRC",
          "PSTR")

data[, qual] <- lapply(data[, qual], as.factor)

# Get the Gower's distance matrix
dist_matrix <- daisy(x = data[, c(qual, quant)],
                     metric = "gower")


data <- cbind(genotypes = rownames(data), data)
row.names(data) <- NULL

# Prepare inputs
counts <- c(I = 31, II = 31, III = 18, IV = 35, V = 40, VI = 17)

mand_accns <-
  c("TMe-34", "TMe-3423", "TMe-2018", "TMe-801", "TMe-551")

gp_vec <- setNames(as.character(data[, "Cluster"]), data[, "genotypes"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by centrality based methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Medoid-like Representative Sampling by Minimal Mean Distance
sel_mean_medoid_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "mean.medoid")
sel_mean_medoid_out
#> $I
#>  [1] "TMe-1696" "TMe-2251" "TMe-3003" "TMe-3202" "TMe-2388" "TMe-1226"
#>  [7] "TMe-1823" "TMe-3416" "TMe-2255" "TMe-990"  "TMe-677"  "TMe-2152"
#> [13] "TMe-1717" "TMe-1940" "TMe-2083" "TMe-1914" "TMe-898"  "TMe-2066"
#> [19] "TMe-3145" "TMe-1448" "TMe-3490" "TMe-1344" "TMe-486"  "TMe-2949"
#> [25] "TMe-3127" "TMe-716"  "TMe-933"  "TMe-2964" "TMe-3001" "TMe-3477"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-433"  "TMe-1969" "TMe-3203" "TMe-3365" "TMe-2242"
#>  [7] "TMe-3617" "TMe-1833" "TMe-377"  "TMe-3540" "TMe-455"  "TMe-2715"
#> [13] "TMe-1323" "TMe-2935" "TMe-1187" "TMe-3671" "TMe-3239" "TMe-2950"
#> [19] "TMe-3146" "TMe-1251" "TMe-3187" "TMe-2765" "TMe-1409" "TMe-2951"
#> [25] "TMe-1766" "TMe-1271" "TMe-1795" "TMe-2705" "TMe-1172" "TMe-2814"
#> [31] "TMe-3642"
#> 
#> $III
#>  [1] "TMe-2154" "TMe-2285" "TMe-3113" "TMe-2169" "TMe-3069" "TMe-3118"
#>  [7] "TMe-304"  "TMe-1138" "TMe-1725" "TMe-3326" "TMe-1836" "TMe-237" 
#> [13] "TMe-1797" "TMe-617"  "TMe-3324" "TMe-3028" "TMe-1675" "TMe-1792"
#> 
#> $IV
#>  [1] "TMe-1726" "TMe-2340" "TMe-3231" "TMe-1652" "TMe-3161" "TMe-2332"
#>  [7] "TMe-3327" "TMe-1867" "TMe-1211" "TMe-2363" "TMe-3763" "TMe-403" 
#> [13] "TMe-3233" "TMe-1828" "TMe-3562" "TMe-1167" "TMe-1651" "TMe-1553"
#> [19] "TMe-3218" "TMe-3168" "TMe-2337" "TMe-3412" "TMe-1118" "TMe-3212"
#> [25] "TMe-2210" "TMe-1948" "TMe-3214" "TMe-3494" "TMe-1996" "TMe-3660"
#> [31] "TMe-2390" "TMe-3232" "TMe-380"  "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-456"  "TMe-474"  "TMe-1414" "TMe-574"  "TMe-2400"
#>  [7] "TMe-2845" "TMe-340"  "TMe-363"  "TMe-1268" "TMe-1160" "TMe-945" 
#> [13] "TMe-565"  "TMe-1234" "TMe-312"  "TMe-1979" "TMe-417"  "TMe-473" 
#> [19] "TMe-359"  "TMe-472"  "TMe-685"  "TMe-647"  "TMe-1862" "TMe-211" 
#> [25] "TMe-424"  "TMe-2271" "TMe-1183" "TMe-1622" "TMe-479"  "TMe-645" 
#> [31] "TMe-245"  "TMe-168"  "TMe-430"  "TMe-1559" "TMe-939"  "TMe-1760"
#> [37] "TMe-747"  "TMe-1382" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-662"  "TMe-791"  "TMe-1847" "TMe-1164" "TMe-907" 
#>  [7] "TMe-1704" "TMe-222"  "TMe-1256" "TMe-985"  "TMe-726"  "TMe-1548"
#> [13] "TMe-1074" "TMe-1481" "TMe-809"  "TMe-720"  "TMe-752" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_mean_medoid_out,
          use.names = FALSE)) +
  labs(title = "mean.medoid")


# Medoid-like Representative Sampling by Minimal Median Distance
sel_median_medoid_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "median.medoid")
sel_median_medoid_out
#> $I
#>  [1] "TMe-3202" "TMe-1696" "TMe-2251" "TMe-3003" "TMe-677"  "TMe-1344"
#>  [7] "TMe-1823" "TMe-2255" "TMe-1226" "TMe-2388" "TMe-3416" "TMe-1451"
#> [13] "TMe-898"  "TMe-3127" "TMe-486"  "TMe-990"  "TMe-1940" "TMe-3145"
#> [19] "TMe-2066" "TMe-2949" "TMe-3001" "TMe-1914" "TMe-1717" "TMe-1955"
#> [25] "TMe-910"  "TMe-2083" "TMe-1448" "TMe-2152" "TMe-716"  "TMe-3142"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3365" "TMe-3495" "TMe-3203" "TMe-3540" "TMe-1833" "TMe-2935"
#>  [7] "TMe-1795" "TMe-2715" "TMe-377"  "TMe-455"  "TMe-433"  "TMe-1766"
#> [13] "TMe-2705" "TMe-2242" "TMe-1271" "TMe-1969" "TMe-2814" "TMe-2765"
#> [19] "TMe-3671" "TMe-1323" "TMe-1187" "TMe-3187" "TMe-2950" "TMe-1172"
#> [25] "TMe-1409" "TMe-2611" "TMe-3239" "TMe-2951" "TMe-2414" "TMe-1505"
#> [31] "TMe-1251"
#> 
#> $III
#>  [1] "TMe-2285" "TMe-3118" "TMe-2154" "TMe-3069" "TMe-1836" "TMe-3113"
#>  [7] "TMe-2169" "TMe-3028" "TMe-1138" "TMe-304"  "TMe-1725" "TMe-3317"
#> [13] "TMe-3324" "TMe-237"  "TMe-1797" "TMe-1421" "TMe-3731" "TMe-3085"
#> 
#> $IV
#>  [1] "TMe-3161" "TMe-3231" "TMe-1726" "TMe-1652" "TMe-2340" "TMe-3562"
#>  [7] "TMe-3212" "TMe-1211" "TMe-1828" "TMe-3327" "TMe-2363" "TMe-1167"
#> [13] "TMe-3660" "TMe-1867" "TMe-3494" "TMe-2337" "TMe-1996" "TMe-1948"
#> [19] "TMe-3412" "TMe-3763" "TMe-2332" "TMe-403"  "TMe-3233" "TMe-3602"
#> [25] "TMe-1553" "TMe-380"  "TMe-734"  "TMe-3218" "TMe-2210" "TMe-1651"
#> [31] "TMe-1118" "TMe-1765" "TMe-1325" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-456"  "TMe-474"  "TMe-574"  "TMe-2845" "TMe-1414"
#>  [7] "TMe-2400" "TMe-565"  "TMe-340"  "TMe-1268" "TMe-312"  "TMe-417" 
#> [13] "TMe-472"  "TMe-363"  "TMe-1160" "TMe-359"  "TMe-945"  "TMe-1979"
#> [19] "TMe-1234" "TMe-647"  "TMe-685"  "TMe-211"  "TMe-747"  "TMe-473" 
#> [25] "TMe-2271" "TMe-745"  "TMe-645"  "TMe-1183" "TMe-1622" "TMe-1760"
#> [31] "TMe-1559" "TMe-1862" "TMe-1308" "TMe-245"  "TMe-479"  "TMe-2361"
#> [37] "TMe-1273" "TMe-1715" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-907"  "TMe-1847" "TMe-1164" "TMe-662"  "TMe-1704"
#>  [7] "TMe-791"  "TMe-1256" "TMe-809"  "TMe-1110" "TMe-222"  "TMe-985" 
#> [13] "TMe-1548" "TMe-1945" "TMe-726"  "TMe-1074" "TMe-683" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_median_medoid_out,
          use.names = FALSE)) +
  labs(title = "median.medoid")


# Representative Sampling by Proximity to Group Centroid
sel_group_centroid_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "nearest.centroid")
sel_group_centroid_out
#> $I
#>  [1] "TMe-1696" "TMe-2251" "TMe-2964" "TMe-3202" "TMe-2388" "TMe-1226"
#>  [7] "TMe-1823" "TMe-3416" "TMe-2255" "TMe-677"  "TMe-990"  "TMe-2152"
#> [13] "TMe-1717" "TMe-2083" "TMe-1940" "TMe-1914" "TMe-898"  "TMe-2066"
#> [19] "TMe-3145" "TMe-1448" "TMe-3490" "TMe-486"  "TMe-1344" "TMe-2949"
#> [25] "TMe-3127" "TMe-2943" "TMe-716"  "TMe-933"  "TMe-3001" "TMe-3477"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-433"  "TMe-1969" "TMe-3203" "TMe-3365" "TMe-2242"
#>  [7] "TMe-3617" "TMe-1833" "TMe-377"  "TMe-3540" "TMe-455"  "TMe-2715"
#> [13] "TMe-1323" "TMe-2935" "TMe-1187" "TMe-3671" "TMe-3239" "TMe-2950"
#> [19] "TMe-3146" "TMe-1251" "TMe-3187" "TMe-2765" "TMe-1409" "TMe-2951"
#> [25] "TMe-1766" "TMe-1271" "TMe-1795" "TMe-2705" "TMe-1172" "TMe-2814"
#> [31] "TMe-3642"
#> 
#> $III
#>  [1] "TMe-2154" "TMe-2285" "TMe-3113" "TMe-2169" "TMe-3069" "TMe-3118"
#>  [7] "TMe-304"  "TMe-1138" "TMe-1725" "TMe-3326" "TMe-1836" "TMe-237" 
#> [13] "TMe-1797" "TMe-617"  "TMe-3324" "TMe-3028" "TMe-1675" "TMe-1792"
#> 
#> $IV
#>  [1] "TMe-1700" "TMe-2332" "TMe-1597" "TMe-3144" "TMe-3218" "TMe-2036"
#>  [7] "TMe-3297" "TMe-1820" "TMe-1179" "TMe-2247" "TMe-3760" "TMe-403" 
#> [13] "TMe-3231" "TMe-3541" "TMe-1652" "TMe-1166" "TMe-3206" "TMe-1580"
#> [19] "TMe-1526" "TMe-3161" "TMe-2318" "TMe-1106" "TMe-2172" "TMe-3390"
#> [25] "TMe-3198" "TMe-1923" "TMe-3451" "TMe-3044" "TMe-3225" "TMe-3619"
#> [31] "TMe-1987" "TMe-1651" "TMe-380"  "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-456"  "TMe-474"  "TMe-1411" "TMe-565"  "TMe-2355"
#>  [7] "TMe-2769" "TMe-340"  "TMe-363"  "TMe-1257" "TMe-1127" "TMe-939" 
#> [13] "TMe-551"  "TMe-1220" "TMe-287"  "TMe-1963" "TMe-417"  "TMe-473" 
#> [19] "TMe-472"  "TMe-359"  "TMe-669"  "TMe-197"  "TMe-1859" "TMe-623" 
#> [25] "TMe-2003" "TMe-424"  "TMe-1534" "TMe-1160" "TMe-479"  "TMe-217" 
#> [31] "TMe-629"  "TMe-142"  "TMe-1551" "TMe-370"  "TMe-932"  "TMe-1730"
#> [37] "TMe-731"  "TMe-2339" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-662"  "TMe-791"  "TMe-1847" "TMe-1164" "TMe-907" 
#>  [7] "TMe-1704" "TMe-222"  "TMe-1256" "TMe-985"  "TMe-726"  "TMe-1548"
#> [13] "TMe-1074" "TMe-1481" "TMe-809"  "TMe-720"  "TMe-752" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_group_centroid_out,
          use.names = FALSE)) +
  labs(title = "nearest.centroid")


# Representative Sampling by Proximity to Group Spatial Median
sel_group_median_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "nearest.median")
sel_group_median_out
#> $I
#>  [1] "TMe-1696" "TMe-2251" "TMe-2964" "TMe-3202" "TMe-2388" "TMe-1226"
#>  [7] "TMe-1823" "TMe-3416" "TMe-2255" "TMe-2152" "TMe-677"  "TMe-1717"
#> [13] "TMe-990"  "TMe-2083" "TMe-1940" "TMe-2066" "TMe-898"  "TMe-3145"
#> [19] "TMe-1914" "TMe-1448" "TMe-486"  "TMe-1344" "TMe-3127" "TMe-933" 
#> [25] "TMe-3490" "TMe-716"  "TMe-2949" "TMe-3001" "TMe-2943" "TMe-1955"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-433"  "TMe-1969" "TMe-3203" "TMe-2242" "TMe-3365"
#>  [7] "TMe-1833" "TMe-3540" "TMe-3617" "TMe-377"  "TMe-455"  "TMe-2715"
#> [13] "TMe-1323" "TMe-2935" "TMe-3671" "TMe-1187" "TMe-3239" "TMe-3146"
#> [19] "TMe-1251" "TMe-2950" "TMe-2765" "TMe-1409" "TMe-3187" "TMe-1766"
#> [25] "TMe-1795" "TMe-2951" "TMe-2705" "TMe-1271" "TMe-2814" "TMe-1172"
#> [31] "TMe-1461"
#> 
#> $III
#>  [1] "TMe-2154" "TMe-3113" "TMe-2285" "TMe-2169" "TMe-3069" "TMe-3118"
#>  [7] "TMe-1725" "TMe-304"  "TMe-1138" "TMe-3326" "TMe-1836" "TMe-237" 
#> [13] "TMe-1797" "TMe-3324" "TMe-617"  "TMe-3028" "TMe-1792" "TMe-1675"
#> 
#> $IV
#>  [1] "TMe-1700" "TMe-2332" "TMe-1597" "TMe-3144" "TMe-3218" "TMe-2036"
#>  [7] "TMe-3297" "TMe-1820" "TMe-1179" "TMe-2247" "TMe-403"  "TMe-3760"
#> [13] "TMe-3231" "TMe-3541" "TMe-1652" "TMe-1580" "TMe-3206" "TMe-1526"
#> [19] "TMe-1166" "TMe-2318" "TMe-1106" "TMe-2172" "TMe-3390" "TMe-1923"
#> [25] "TMe-3198" "TMe-3161" "TMe-1651" "TMe-3451" "TMe-3225" "TMe-380" 
#> [31] "TMe-3619" "TMe-3044" "TMe-1987" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-456"  "TMe-474"  "TMe-1411" "TMe-565"  "TMe-2769"
#>  [7] "TMe-2355" "TMe-340"  "TMe-363"  "TMe-1127" "TMe-1257" "TMe-939" 
#> [13] "TMe-551"  "TMe-287"  "TMe-1220" "TMe-1963" "TMe-473"  "TMe-417" 
#> [19] "TMe-359"  "TMe-472"  "TMe-197"  "TMe-669"  "TMe-1859" "TMe-1534"
#> [25] "TMe-623"  "TMe-1160" "TMe-2003" "TMe-424"  "TMe-629"  "TMe-217" 
#> [31] "TMe-479"  "TMe-2339" "TMe-370"  "TMe-932"  "TMe-1551" "TMe-731" 
#> [37] "TMe-764"  "TMe-1694" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-662"  "TMe-791"  "TMe-1847" "TMe-1164" "TMe-907" 
#>  [7] "TMe-1704" "TMe-222"  "TMe-1256" "TMe-985"  "TMe-726"  "TMe-1074"
#> [13] "TMe-1548" "TMe-1481" "TMe-809"  "TMe-720"  "TMe-1110"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_group_median_out,
          use.names = FALSE)) +
  labs(title = "nearest.median")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by peripheral/extremity based methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Peripheral Sampling by Maximal Mean Distance
sel_mean_peripheral_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "mean.peripheral")
sel_mean_peripheral_out
#> $I
#>  [1] "TMe-2967" "TMe-3685" "TMe-3705" "TMe-3292" "TMe-3223" "TMe-3667"
#>  [7] "TMe-41"   "TMe-3475" "TMe-606"  "TMe-1425" "TMe-3282" "TMe-3319"
#> [13] "TMe-3392" "TMe-2965" "TMe-3437" "TMe-1564" "TMe-2943" "TMe-2955"
#> [19] "TMe-2993" "TMe-756"  "TMe-2069" "TMe-717"  "TMe-3415" "TMe-3396"
#> [25] "TMe-3389" "TMe-500"  "TMe-3065" "TMe-2996" "TMe-867"  "TMe-3641"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3466" "TMe-3338" "TMe-47"   "TMe-2860" "TMe-1919" "TMe-3605"
#>  [7] "TMe-1907" "TMe-86"   "TMe-1107" "TMe-2128" "TMe-369"  "TMe-2033"
#> [13] "TMe-2952" "TMe-2056" "TMe-3141" "TMe-3805" "TMe-2568" "TMe-2329"
#> [19] "TMe-1385" "TMe-1137" "TMe-3766" "TMe-3800" "TMe-674"  "TMe-3210"
#> [25] "TMe-3200" "TMe-2352" "TMe-2891" "TMe-3690" "TMe-3140" "TMe-509" 
#> [31] "TMe-1619"
#> 
#> $III
#>  [1] "TMe-3596" "TMe-2926" "TMe-1897" "TMe-3274" "TMe-3565" "TMe-267" 
#>  [7] "TMe-3804" "TMe-3299" "TMe-785"  "TMe-234"  "TMe-635"  "TMe-261" 
#> [13] "TMe-1868" "TMe-3445" "TMe-35"   "TMe-3663" "TMe-116"  "TMe-773" 
#> 
#> $IV
#>  [1] "TMe-2809" "TMe-241"  "TMe-1434" "TMe-550"  "TMe-812"  "TMe-761" 
#>  [7] "TMe-3417" "TMe-2924" "TMe-3442" "TMe-1988" "TMe-1525" "TMe-698" 
#> [13] "TMe-1020" "TMe-3730" "TMe-1078" "TMe-3054" "TMe-608"  "TMe-3273"
#> [19] "TMe-3573" "TMe-2043" "TMe-226"  "TMe-180"  "TMe-2956" "TMe-3275"
#> [25] "TMe-2567" "TMe-154"  "TMe-3545" "TMe-421"  "TMe-1229" "TMe-2971"
#> [31] "TMe-317"  "TMe-3025" "TMe-3707" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-3034" "TMe-2121" "TMe-3628" "TMe-419"  "TMe-3736" "TMe-3185"
#>  [7] "TMe-536"  "TMe-2853" "TMe-378"  "TMe-712"  "TMe-2355" "TMe-2589"
#> [13] "TMe-2050" "TMe-1458" "TMe-2531" "TMe-616"  "TMe-723"  "TMe-3356"
#> [19] "TMe-1381" "TMe-1099" "TMe-1455" "TMe-3332" "TMe-399"  "TMe-700" 
#> [25] "TMe-782"  "TMe-2612" "TMe-603"  "TMe-277"  "TMe-861"  "TMe-2953"
#> [31] "TMe-284"  "TMe-1391" "TMe-2907" "TMe-2820" "TMe-2296" "TMe-3411"
#> [37] "TMe-2413" "TMe-3408" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-3095" "TMe-3549" "TMe-2983" "TMe-2035" "TMe-3116" "TMe-1416"
#>  [7] "TMe-693"  "TMe-2963" "TMe-1232" "TMe-1769" "TMe-1403" "TMe-2064"
#> [13] "TMe-1875" "TMe-1124" "TMe-1383" "TMe-1518" "TMe-1428"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_mean_peripheral_out,
          use.names = FALSE)) +
  labs(title = "mean.peripheral")


# Peripheral Sampling by Maximal Median Distance
sel_median_peripheral_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "median.peripheral")
sel_median_peripheral_out
#> $I
#>  [1] "TMe-3685" "TMe-2967" "TMe-3705" "TMe-3223" "TMe-3292" "TMe-3667"
#>  [7] "TMe-3475" "TMe-41"   "TMe-3282" "TMe-3319" "TMe-1425" "TMe-606" 
#> [13] "TMe-3392" "TMe-2993" "TMe-2955" "TMe-2943" "TMe-2965" "TMe-1564"
#> [19] "TMe-756"  "TMe-3437" "TMe-2069" "TMe-3641" "TMe-3389" "TMe-3415"
#> [25] "TMe-3065" "TMe-3396" "TMe-717"  "TMe-2779" "TMe-3249" "TMe-3266"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3466" "TMe-47"   "TMe-3338" "TMe-1919" "TMe-2860" "TMe-3605"
#>  [7] "TMe-1907" "TMe-1107" "TMe-2128" "TMe-86"   "TMe-369"  "TMe-2056"
#> [13] "TMe-2952" "TMe-2033" "TMe-674"  "TMe-3141" "TMe-1385" "TMe-1137"
#> [19] "TMe-3805" "TMe-2568" "TMe-3690" "TMe-3800" "TMe-1619" "TMe-2329"
#> [25] "TMe-2891" "TMe-3766" "TMe-3210" "TMe-2352" "TMe-509"  "TMe-3200"
#> [31] "TMe-1242"
#> 
#> $III
#>  [1] "TMe-3596" "TMe-2926" "TMe-1897" "TMe-3274" "TMe-3565" "TMe-3299"
#>  [7] "TMe-267"  "TMe-116"  "TMe-234"  "TMe-261"  "TMe-785"  "TMe-3663"
#> [13] "TMe-3804" "TMe-773"  "TMe-3631" "TMe-35"   "TMe-635"  "TMe-1868"
#> 
#> $IV
#>  [1] "TMe-2809" "TMe-241"  "TMe-1434" "TMe-550"  "TMe-812"  "TMe-761" 
#>  [7] "TMe-2924" "TMe-3417" "TMe-1988" "TMe-1078" "TMe-3442" "TMe-1020"
#> [13] "TMe-608"  "TMe-698"  "TMe-1525" "TMe-3273" "TMe-3730" "TMe-3573"
#> [19] "TMe-3054" "TMe-317"  "TMe-226"  "TMe-2043" "TMe-2971" "TMe-2956"
#> [25] "TMe-2567" "TMe-1419" "TMe-1229" "TMe-154"  "TMe-180"  "TMe-3545"
#> [31] "TMe-650"  "TMe-421"  "TMe-3707" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-2121" "TMe-3034" "TMe-3628" "TMe-419"  "TMe-3185" "TMe-3736"
#>  [7] "TMe-536"  "TMe-1458" "TMe-378"  "TMe-2853" "TMe-1455" "TMe-1381"
#> [13] "TMe-2355" "TMe-2589" "TMe-712"  "TMe-2050" "TMe-723"  "TMe-616" 
#> [19] "TMe-3356" "TMe-3332" "TMe-2531" "TMe-1099" "TMe-700"  "TMe-1391"
#> [25] "TMe-2612" "TMe-782"  "TMe-277"  "TMe-399"  "TMe-2296" "TMe-2953"
#> [31] "TMe-2907" "TMe-603"  "TMe-2820" "TMe-861"  "TMe-284"  "TMe-2790"
#> [37] "TMe-832"  "TMe-2413" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-3095" "TMe-3549" "TMe-1403" "TMe-2983" "TMe-3116" "TMe-693" 
#>  [7] "TMe-1416" "TMe-2035" "TMe-1232" "TMe-1769" "TMe-2963" "TMe-2064"
#> [13] "TMe-1875" "TMe-1124" "TMe-1383" "TMe-1518" "TMe-1775"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_median_peripheral_out,
          use.names = FALSE)) +
  labs(title = "median.peripheral")


# Peripheral Sampling by Maximal Eccentricity
sel_eccentricity_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "eccentricity")
sel_eccentricity_out
#> $I
#>  [1] "TMe-606"  "TMe-3667" "TMe-3030" "TMe-3705" "TMe-3685" "TMe-1533"
#>  [7] "TMe-3480" "TMe-3437" "TMe-2943" "TMe-132"  "TMe-3719" "TMe-1564"
#> [13] "TMe-3292" "TMe-41"   "TMe-2914" "TMe-3698" "TMe-2967" "TMe-3169"
#> [19] "TMe-3208" "TMe-2383" "TMe-1347" "TMe-2004" "TMe-815"  "TMe-3223"
#> [25] "TMe-2216" "TMe-896"  "TMe-1621" "TMe-3110" "TMe-2965" "TMe-1468"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-369"  "TMe-3466" "TMe-2056" "TMe-3605" "TMe-86"   "TMe-1107"
#>  [7] "TMe-47"   "TMe-774"  "TMe-160"  "TMe-674"  "TMe-404"  "TMe-3141"
#> [13] "TMe-3766" "TMe-1242" "TMe-196"  "TMe-2329" "TMe-289"  "TMe-1385"
#> [19] "TMe-2204" "TMe-1907" "TMe-3771" "TMe-176"  "TMe-1739" "TMe-3338"
#> [25] "TMe-478"  "TMe-509"  "TMe-3406" "TMe-3009" "TMe-2860" "TMe-2352"
#> [31] "TMe-171" 
#> 
#> $III
#>  [1] "TMe-878"  "TMe-3596" "TMe-785"  "TMe-2926" "TMe-3128" "TMe-1910"
#>  [7] "TMe-1897" "TMe-635"  "TMe-1863" "TMe-267"  "TMe-200"  "TMe-13"  
#> [13] "TMe-1443" "TMe-1804" "TMe-773"  "TMe-1889" "TMe-3274" "TMe-1790"
#> 
#> $IV
#>  [1] "TMe-2809" "TMe-609"  "TMe-1335" "TMe-241"  "TMe-596"  "TMe-1179"
#>  [7] "TMe-280"  "TMe-1975" "TMe-421"  "TMe-3417" "TMe-1525" "TMe-3386"
#> [13] "TMe-1988" "TMe-2788" "TMe-3040" "TMe-1095" "TMe-427"  "TMe-2833"
#> [19] "TMe-2956" "TMe-1434" "TMe-3025" "TMe-802"  "TMe-226"  "TMe-761" 
#> [25] "TMe-2567" "TMe-3273" "TMe-1336" "TMe-608"  "TMe-812"  "TMe-3316"
#> [31] "TMe-1078" "TMe-3581" "TMe-3576" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-399"  "TMe-3034" "TMe-2531" "TMe-277"  "TMe-1534" "TMe-700" 
#>  [7] "TMe-2121" "TMe-294"  "TMe-3628" "TMe-378"  "TMe-2853" "TMe-419" 
#> [13] "TMe-603"  "TMe-2933" "TMe-2050" "TMe-2907" "TMe-627"  "TMe-2296"
#> [19] "TMe-832"  "TMe-3411" "TMe-1311" "TMe-1399" "TMe-2790" "TMe-1011"
#> [25] "TMe-2589" "TMe-1307" "TMe-2003" "TMe-1294" "TMe-731"  "TMe-1158"
#> [31] "TMe-336"  "TMe-1099" "TMe-1455" "TMe-1257" "TMe-1300" "TMe-2435"
#> [37] "TMe-2953" "TMe-1730" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1428" "TMe-2963" "TMe-2035" "TMe-1232" "TMe-1403" "TMe-2983"
#>  [7] "TMe-1079" "TMe-3549" "TMe-1566" "TMe-1900" "TMe-3095" "TMe-705" 
#> [13] "TMe-1646" "TMe-1858" "TMe-836"  "TMe-3387" "TMe-690" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_eccentricity_out,
          use.names = FALSE)) +
  labs(title = "eccentricity")


# Peripheral Sampling by Maximal Farness Centrality
sel_far_cent_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "farness.centrality")
sel_far_cent_out
#> $I
#>  [1] "TMe-2967" "TMe-3685" "TMe-3705" "TMe-3292" "TMe-3223" "TMe-3667"
#>  [7] "TMe-41"   "TMe-3475" "TMe-606"  "TMe-1425" "TMe-3282" "TMe-3319"
#> [13] "TMe-3392" "TMe-2965" "TMe-3437" "TMe-1564" "TMe-2943" "TMe-2955"
#> [19] "TMe-2993" "TMe-756"  "TMe-2069" "TMe-717"  "TMe-3415" "TMe-3396"
#> [25] "TMe-3389" "TMe-500"  "TMe-3065" "TMe-2996" "TMe-867"  "TMe-3641"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3466" "TMe-3338" "TMe-47"   "TMe-2860" "TMe-1919" "TMe-3605"
#>  [7] "TMe-1907" "TMe-86"   "TMe-1107" "TMe-2128" "TMe-369"  "TMe-2033"
#> [13] "TMe-2952" "TMe-2056" "TMe-3141" "TMe-3805" "TMe-2568" "TMe-2329"
#> [19] "TMe-1385" "TMe-1137" "TMe-3766" "TMe-3800" "TMe-674"  "TMe-3210"
#> [25] "TMe-3200" "TMe-2352" "TMe-2891" "TMe-3690" "TMe-3140" "TMe-509" 
#> [31] "TMe-1619"
#> 
#> $III
#>  [1] "TMe-3596" "TMe-2926" "TMe-1897" "TMe-3274" "TMe-3565" "TMe-267" 
#>  [7] "TMe-3804" "TMe-3299" "TMe-785"  "TMe-234"  "TMe-635"  "TMe-261" 
#> [13] "TMe-1868" "TMe-3445" "TMe-35"   "TMe-3663" "TMe-116"  "TMe-773" 
#> 
#> $IV
#>  [1] "TMe-2809" "TMe-241"  "TMe-1434" "TMe-550"  "TMe-812"  "TMe-761" 
#>  [7] "TMe-3417" "TMe-2924" "TMe-3442" "TMe-1988" "TMe-1525" "TMe-698" 
#> [13] "TMe-1020" "TMe-3730" "TMe-1078" "TMe-3054" "TMe-608"  "TMe-3273"
#> [19] "TMe-3573" "TMe-2043" "TMe-226"  "TMe-180"  "TMe-2956" "TMe-3275"
#> [25] "TMe-2567" "TMe-154"  "TMe-3545" "TMe-421"  "TMe-1229" "TMe-2971"
#> [31] "TMe-317"  "TMe-3025" "TMe-3707" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-3034" "TMe-2121" "TMe-3628" "TMe-419"  "TMe-3736" "TMe-3185"
#>  [7] "TMe-536"  "TMe-2853" "TMe-378"  "TMe-712"  "TMe-2355" "TMe-2589"
#> [13] "TMe-2050" "TMe-1458" "TMe-2531" "TMe-616"  "TMe-723"  "TMe-3356"
#> [19] "TMe-1381" "TMe-1099" "TMe-1455" "TMe-3332" "TMe-399"  "TMe-700" 
#> [25] "TMe-782"  "TMe-2612" "TMe-603"  "TMe-277"  "TMe-861"  "TMe-2953"
#> [31] "TMe-284"  "TMe-1391" "TMe-2907" "TMe-2820" "TMe-2296" "TMe-3411"
#> [37] "TMe-2413" "TMe-3408" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-3095" "TMe-3549" "TMe-2983" "TMe-2035" "TMe-3116" "TMe-1416"
#>  [7] "TMe-693"  "TMe-2963" "TMe-1232" "TMe-1769" "TMe-1403" "TMe-2064"
#> [13] "TMe-1875" "TMe-1124" "TMe-1383" "TMe-1518" "TMe-1428"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_far_cent_out,
          use.names = FALSE)) +
  labs(title = "farness.centrality")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by space-Filling/coverage methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Space-Filling Sampling via the Kennard-Stone Algorithm
sel_ks_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "kennard.stone")
sel_ks_out
#> $I
#>  [1] "TMe-3268" "TMe-937"  "TMe-3623" "TMe-1621" "TMe-3481" "TMe-3719"
#>  [7] "TMe-501"  "TMe-2984" "TMe-3163" "TMe-3208" "TMe-469"  "TMe-815" 
#> [13] "TMe-1568" "TMe-3330" "TMe-3521" "TMe-438"  "TMe-1347" "TMe-3145"
#> [19] "TMe-3219" "TMe-2462" "TMe-3132" "TMe-3292" "TMe-3392" "TMe-3518"
#> [25] "TMe-33"   "TMe-446"  "TMe-736"  "TMe-778"  "TMe-867"  "TMe-1096"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3571" "TMe-589"  "TMe-890"  "TMe-660"  "TMe-1409" "TMe-706" 
#>  [7] "TMe-539"  "TMe-1831" "TMe-3447" "TMe-3467" "TMe-40"   "TMe-1619"
#> [13] "TMe-796"  "TMe-2200" "TMe-2851" "TMe-2978" "TMe-3009" "TMe-3222"
#> [19] "TMe-774"  "TMe-2"    "TMe-85"   "TMe-1051" "TMe-1754" "TMe-1864"
#> [25] "TMe-2304" "TMe-2797" "TMe-3114" "TMe-3557" "TMe-196"  "TMe-1459"
#> [31] "TMe-1616"
#> 
#> $III
#>  [1] "TMe-3715" "TMe-1421" "TMe-3287" "TMe-1593" "TMe-14"   "TMe-434" 
#>  [7] "TMe-141"  "TMe-70"   "TMe-2926" "TMe-3008" "TMe-773"  "TMe-2270"
#> [13] "TMe-203"  "TMe-381"  "TMe-2086" "TMe-2502" "TMe-1230" "TMe-425" 
#> 
#> $IV
#>  [1] "TMe-3340" "TMe-90"   "TMe-1776" "TMe-519"  "TMe-27"   "TMe-12"  
#>  [7] "TMe-2458" "TMe-1526" "TMe-3144" "TMe-3253" "TMe-107"  "TMe-526" 
#> [13] "TMe-2954" "TMe-1419" "TMe-3276" "TMe-57"   "TMe-204"  "TMe-1095"
#> [19] "TMe-1155" "TMe-2764" "TMe-3191" "TMe-81"   "TMe-887"  "TMe-3412"
#> [25] "TMe-3639" "TMe-15"   "TMe-3484" "TMe-154"  "TMe-396"  "TMe-760" 
#> [31] "TMe-1246" "TMe-1814" "TMe-2039" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-1551" "TMe-373"  "TMe-600"  "TMe-1345" "TMe-1733" "TMe-1411"
#>  [7] "TMe-487"  "TMe-1534" "TMe-886"  "TMe-1160" "TMe-1418" "TMe-2213"
#> [13] "TMe-167"  "TMe-217"  "TMe-655"  "TMe-850"  "TMe-399"  "TMe-997" 
#> [19] "TMe-1500" "TMe-419"  "TMe-456"  "TMe-1268" "TMe-1373" "TMe-1880"
#> [25] "TMe-2121" "TMe-1322" "TMe-347"  "TMe-378"  "TMe-2530" "TMe-168" 
#> [31] "TMe-712"  "TMe-1788" "TMe-3311" "TMe-423"  "TMe-474"  "TMe-603" 
#> [37] "TMe-1037" "TMe-1963" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1509" "TMe-1985" "TMe-1481" "TMe-693"  "TMe-1233" "TMe-696" 
#>  [7] "TMe-2791" "TMe-222"  "TMe-3549" "TMe-631"  "TMe-1483" "TMe-531" 
#> [13] "TMe-985"  "TMe-1261" "TMe-1413" "TMe-1428" "TMe-662" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_ks_out,
          use.names = FALSE)) +
  labs(title = "kennard.stone")


# Space-Filling Sampling via the DUPLEX Algorithm
sel_duplex_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "duplex")
sel_duplex_out
#> $I
#>  [1] "TMe-3268" "TMe-937"  "TMe-3623" "TMe-1621" "TMe-3481" "TMe-3719"
#>  [7] "TMe-501"  "TMe-2984" "TMe-3163" "TMe-3208" "TMe-469"  "TMe-1568"
#> [13] "TMe-3330" "TMe-438"  "TMe-1347" "TMe-3145" "TMe-2462" "TMe-2909"
#> [19] "TMe-3132" "TMe-3292" "TMe-3518" "TMe-33"   "TMe-446"  "TMe-736" 
#> [25] "TMe-778"  "TMe-867"  "TMe-1096" "TMe-1448" "TMe-1823" "TMe-2388"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3571" "TMe-589"  "TMe-890"  "TMe-660"  "TMe-1409" "TMe-706" 
#>  [7] "TMe-539"  "TMe-2715" "TMe-3447" "TMe-40"   "TMe-1619" "TMe-796" 
#> [13] "TMe-1137" "TMe-2851" "TMe-85"   "TMe-339"  "TMe-377"  "TMe-1051"
#> [19] "TMe-1250" "TMe-1732" "TMe-2304" "TMe-2797" "TMe-1459" "TMe-3443"
#> [25] "TMe-2"    "TMe-369"  "TMe-1005" "TMe-1754" "TMe-2054" "TMe-2862"
#> [31] "TMe-53"  
#> 
#> $III
#>  [1] "TMe-3715" "TMe-1421" "TMe-3287" "TMe-1593" "TMe-14"   "TMe-2926"
#>  [7] "TMe-141"  "TMe-2203" "TMe-1910" "TMe-6"    "TMe-206"  "TMe-3362"
#> [13] "TMe-70"   "TMe-3143" "TMe-1176" "TMe-1868" "TMe-2161" "TMe-3147"
#> 
#> $IV
#>  [1] "TMe-3340" "TMe-90"   "TMe-1776" "TMe-519"  "TMe-27"   "TMe-12"  
#>  [7] "TMe-2458" "TMe-1526" "TMe-3144" "TMe-3253" "TMe-526"  "TMe-2954"
#> [13] "TMe-1419" "TMe-1765" "TMe-3276" "TMe-204"  "TMe-834"  "TMe-1155"
#> [19] "TMe-3257" "TMe-887"  "TMe-3639" "TMe-3484" "TMe-396"  "TMe-760" 
#> [25] "TMe-1246" "TMe-3072" "TMe-3435" "TMe-30"   "TMe-207"  "TMe-1651"
#> [31] "TMe-2989" "TMe-3019" "TMe-3494" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-1551" "TMe-373"  "TMe-600"  "TMe-1345" "TMe-1733" "TMe-1411"
#>  [7] "TMe-487"  "TMe-1534" "TMe-1160" "TMe-1418" "TMe-2213" "TMe-217" 
#> [13] "TMe-655"  "TMe-850"  "TMe-399"  "TMe-2015" "TMe-419"  "TMe-1268"
#> [19] "TMe-1373" "TMe-1880" "TMe-748"  "TMe-1322" "TMe-347"  "TMe-2530"
#> [25] "TMe-168"  "TMe-370"  "TMe-712"  "TMe-3311" "TMe-474"  "TMe-603" 
#> [31] "TMe-667"  "TMe-1963" "TMe-771"  "TMe-532"  "TMe-574"  "TMe-945" 
#> [37] "TMe-1381" "TMe-1762" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1509" "TMe-1985" "TMe-1481" "TMe-693"  "TMe-1233" "TMe-696" 
#>  [7] "TMe-2791" "TMe-3549" "TMe-631"  "TMe-1483" "TMe-662"  "TMe-531" 
#> [13] "TMe-985"  "TMe-1261" "TMe-1413" "TMe-1428" "TMe-1678"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_duplex_out,
          use.names = FALSE)) +
  labs(title = "duplex")


# Space-Filling Sampling via the Honigs Algorithm
sel_honigs_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "honigs")
sel_honigs_out
#> $I
#>  [1] "TMe-3667" "TMe-606"  "TMe-2965" "TMe-1423" "TMe-2967" "TMe-3337"
#>  [7] "TMe-3223" "TMe-2513" "TMe-3163" "TMe-3151" "TMe-3478" "TMe-2993"
#> [13] "TMe-2934" "TMe-1786" "TMe-3481" "TMe-1589" "TMe-715"  "TMe-264" 
#> [19] "TMe-1425" "TMe-867"  "TMe-2975" "TMe-3026" "TMe-299"  "TMe-3437"
#> [25] "TMe-579"  "TMe-3471" "TMe-3175" "TMe-132"  "TMe-2462" "TMe-3292"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3466" "TMe-369"  "TMe-3605" "TMe-1919" "TMe-1907" "TMe-3805"
#>  [7] "TMe-47"   "TMe-171"  "TMe-3338" "TMe-289"  "TMe-2329" "TMe-3690"
#> [13] "TMe-3366" "TMe-3141" "TMe-3547" "TMe-431"  "TMe-3150" "TMe-160" 
#> [19] "TMe-3800" "TMe-3209" "TMe-1051" "TMe-2204" "TMe-3530" "TMe-2891"
#> [25] "TMe-2952" "TMe-60"   "TMe-2033" "TMe-3766" "TMe-3222" "TMe-3286"
#> [31] "TMe-509" 
#> 
#> $III
#>  [1] "TMe-3596" "TMe-878"  "TMe-234"  "TMe-425"  "TMe-1897" "TMe-946" 
#>  [7] "TMe-35"   "TMe-3274" "TMe-13"   "TMe-6"    "TMe-1993" "TMe-3299"
#> [13] "TMe-1200" "TMe-3679" "TMe-3715" "TMe-2374" "TMe-3700" "TMe-2086"
#> 
#> $IV
#>  [1] "TMe-609"  "TMe-2809" "TMe-761"  "TMe-3054" "TMe-550"  "TMe-812" 
#>  [7] "TMe-2924" "TMe-225"  "TMe-2971" "TMe-3025" "TMe-3542" "TMe-3273"
#> [13] "TMe-3422" "TMe-656"  "TMe-1020" "TMe-3428" "TMe-3527" "TMe-3275"
#> [19] "TMe-2240" "TMe-3573" "TMe-3707" "TMe-3730" "TMe-317"  "TMe-226" 
#> [25] "TMe-3593" "TMe-2567" "TMe-733"  "TMe-372"  "TMe-241"  "TMe-1511"
#> [31] "TMe-3055" "TMe-1988" "TMe-919"  "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-3034" "TMe-399"  "TMe-712"  "TMe-3185" "TMe-1381" "TMe-3332"
#>  [7] "TMe-277"  "TMe-2853" "TMe-2050" "TMe-832"  "TMe-723"  "TMe-700" 
#> [13] "TMe-2933" "TMe-798"  "TMe-1366" "TMe-3736" "TMe-870"  "TMe-759" 
#> [19] "TMe-487"  "TMe-895"  "TMe-336"  "TMe-1290" "TMe-3628" "TMe-1401"
#> [25] "TMe-208"  "TMe-2413" "TMe-2790" "TMe-1042" "TMe-2905" "TMe-2296"
#> [31] "TMe-2688" "TMe-1026" "TMe-782"  "TMe-2139" "TMe-2973" "TMe-391" 
#> [37] "TMe-536"  "TMe-3356" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-2963" "TMe-1428" "TMe-1875" "TMe-1035" "TMe-1232" "TMe-361" 
#>  [7] "TMe-3095" "TMe-2035" "TMe-1769" "TMe-693"  "TMe-315"  "TMe-1124"
#> [13] "TMe-751"  "TMe-1796" "TMe-531"  "TMe-2791" "TMe-1174"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_honigs_out,
          use.names = FALSE)) +
  labs(title = "honigs")


# Space-Filling Sampling via Farthest-Point (Max-Min) Algorithm
sel_far_pt_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "farthest.point")
sel_far_pt_out
#> $I
#>  [1] "TMe-2216" "TMe-3705" "TMe-2965" "TMe-756"  "TMe-264"  "TMe-3437"
#>  [7] "TMe-1425" "TMe-3266" "TMe-2912" "TMe-1589" "TMe-3415" "TMe-778" 
#> [13] "TMe-2934" "TMe-300"  "TMe-3514" "TMe-1786" "TMe-2943" "TMe-1922"
#> [19] "TMe-3475" "TMe-566"  "TMe-3462" "TMe-3163" "TMe-2993" "TMe-715" 
#> [25] "TMe-410"  "TMe-3478" "TMe-3292" "TMe-1958" "TMe-3625" "TMe-2955"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2352" "TMe-369"  "TMe-3605" "TMe-2056" "TMe-3466" "TMe-1459"
#>  [7] "TMe-2860" "TMe-3690" "TMe-47"   "TMe-3447" "TMe-3210" "TMe-3338"
#> [13] "TMe-2033" "TMe-509"  "TMe-1619" "TMe-768"  "TMe-1919" "TMe-3455"
#> [19] "TMe-3800" "TMe-3406" "TMe-3200" "TMe-1385" "TMe-902"  "TMe-289" 
#> [25] "TMe-3052" "TMe-2952" "TMe-3009" "TMe-648"  "TMe-1907" "TMe-1754"
#> [31] "TMe-2128"
#> 
#> $III
#>  [1] "TMe-2205" "TMe-3565" "TMe-3596" "TMe-635"  "TMe-2977" "TMe-1897"
#>  [7] "TMe-3274" "TMe-2203" "TMe-4"    "TMe-3299" "TMe-3721" "TMe-35"  
#> [13] "TMe-234"  "TMe-123"  "TMe-3445" "TMe-3715" "TMe-3397" "TMe-785" 
#> 
#> $IV
#>  [1] "TMe-3417" "TMe-421"  "TMe-3273" "TMe-698"  "TMe-1700" "TMe-2809"
#>  [7] "TMe-3701" "TMe-2240" "TMe-3054" "TMe-1434" "TMe-2775" "TMe-619" 
#> [13] "TMe-226"  "TMe-3573" "TMe-761"  "TMe-3730" "TMe-824"  "TMe-766" 
#> [19] "TMe-1350" "TMe-3760" "TMe-2839" "TMe-3576" "TMe-383"  "TMe-3297"
#> [25] "TMe-427"  "TMe-184"  "TMe-1278" "TMe-3428" "TMe-967"  "TMe-737" 
#> [31] "TMe-3422" "TMe-2956" "TMe-3055" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-2139" "TMe-700"  "TMe-2531" "TMe-3311" "TMe-3736" "TMe-2853"
#>  [7] "TMe-2863" "TMe-2121" "TMe-194"  "TMe-1500" "TMe-481"  "TMe-655" 
#> [13] "TMe-603"  "TMe-284"  "TMe-336"  "TMe-3185" "TMe-616"  "TMe-1257"
#> [19] "TMe-294"  "TMe-2589" "TMe-3034" "TMe-362"  "TMe-816"  "TMe-582" 
#> [25] "TMe-819"  "TMe-2907" "TMe-3628" "TMe-287"  "TMe-1381" "TMe-1311"
#> [31] "TMe-3411" "TMe-305"  "TMe-378"  "TMe-277"  "TMe-2355" "TMe-1643"
#> [37] "TMe-588"  "TMe-208"  "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1403" "TMe-2035" "TMe-580"  "TMe-693"  "TMe-2534" "TMe-751" 
#>  [7] "TMe-1232" "TMe-1985" "TMe-2791" "TMe-1796" "TMe-2045" "TMe-1079"
#> [13] "TMe-1416" "TMe-310"  "TMe-548"  "TMe-1775" "TMe-1174"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_far_pt_out,
          use.names = FALSE)) +
  labs(title = "farthest.point")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by density based methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Density-Based Sampling by Minimal Nearest-Neighbour Distance
sel_nn_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "nearest.neighbour")
sel_nn_out
#> $I
#>  [1] "TMe-1696" "TMe-1717" "TMe-2976" "TMe-20"   "TMe-2913" "TMe-2916"
#>  [7] "TMe-839"  "TMe-1086" "TMe-1140" "TMe-2084" "TMe-2103" "TMe-2152"
#> [13] "TMe-1423" "TMe-2131" "TMe-438"  "TMe-3424" "TMe-1448" "TMe-2191"
#> [19] "TMe-1221" "TMe-3170" "TMe-3548" "TMe-486"  "TMe-2251" "TMe-990" 
#> [25] "TMe-1823" "TMe-3325" "TMe-3729" "TMe-583"  "TMe-1981" "TMe-3027"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-455"  "TMe-890"  "TMe-3066" "TMe-3571" "TMe-528"  "TMe-930" 
#>  [7] "TMe-1505" "TMe-2814" "TMe-2127" "TMe-664"  "TMe-1251" "TMe-1833"
#> [13] "TMe-3203" "TMe-1409" "TMe-1766" "TMe-2021" "TMe-2412" "TMe-365" 
#> [19] "TMe-2611" "TMe-3533" "TMe-377"  "TMe-2935" "TMe-433"  "TMe-2705"
#> [25] "TMe-2765" "TMe-3187" "TMe-2951" "TMe-1969" "TMe-3642" "TMe-3531"
#> [31] "TMe-3239"
#> 
#> $III
#>  [1] "TMe-2161" "TMe-10"   "TMe-853"  "TMe-1738" "TMe-304"  "TMe-1836"
#>  [7] "TMe-3631" "TMe-3721" "TMe-64"   "TMe-3346" "TMe-1797" "TMe-1910"
#> [13] "TMe-3118" "TMe-617"  "TMe-3113" "TMe-3147" "TMe-3336" "TMe-1790"
#> 
#> $IV
#>  [1] "TMe-2247" "TMe-1139" "TMe-3231" "TMe-3660" "TMe-1118" "TMe-1402"
#>  [7] "TMe-107"  "TMe-72"   "TMe-885"  "TMe-2410" "TMe-57"   "TMe-204" 
#> [13] "TMe-150"  "TMe-1364" "TMe-24"   "TMe-3044" "TMe-3243" "TMe-3779"
#> [19] "TMe-1726" "TMe-1948" "TMe-3357" "TMe-7"    "TMe-380"  "TMe-218" 
#> [25] "TMe-396"  "TMe-3639" "TMe-3494" "TMe-103"  "TMe-3233" "TMe-3240"
#> [31] "TMe-428"  "TMe-1765" "TMe-1374" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-424"  "TMe-430"  "TMe-1559" "TMe-977"  "TMe-1234" "TMe-1572"
#>  [7] "TMe-1879" "TMe-2009" "TMe-2441" "TMe-2361" "TMe-212"  "TMe-464" 
#> [13] "TMe-142"  "TMe-669"  "TMe-1290" "TMe-474"  "TMe-1414" "TMe-618" 
#> [19] "TMe-889"  "TMe-1160" "TMe-1862" "TMe-137"  "TMe-360"  "TMe-1521"
#> [25] "TMe-565"  "TMe-2845" "TMe-168"  "TMe-1188" "TMe-754"  "TMe-330" 
#> [31] "TMe-1367" "TMe-1901" "TMe-1294" "TMe-2439" "TMe-1871" "TMe-1541"
#> [37] "TMe-2213" "TMe-2319" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-720"  "TMe-726"  "TMe-1481" "TMe-2535" "TMe-1074" "TMe-936" 
#>  [7] "TMe-138"  "TMe-2195" "TMe-1945" "TMe-705"  "TMe-985"  "TMe-1164"
#> [13] "TMe-281"  "TMe-1286" "TMe-983"  "TMe-1554" "TMe-2389"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_nn_out,
          use.names = FALSE)) +
  labs(title = "nearest.neighbour")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by cluster based methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Globally Optimal Medoid Sampling via Partitioning Around Medoids (PAM)
sel_pam_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "optim.medoid")
sel_pam_out
#> $I
#>  [1] "TMe-1981" "TMe-3729" "TMe-1777" "TMe-1532" "TMe-3424" "TMe-2103"
#>  [7] "TMe-3465" "TMe-839"  "TMe-2191" "TMe-1341" "TMe-1716" "TMe-3202"
#> [13] "TMe-3003" "TMe-3577" "TMe-756"  "TMe-3425" "TMe-2964" "TMe-882" 
#> [19] "TMe-1823" "TMe-610"  "TMe-2131" "TMe-2251" "TMe-2066" "TMe-2604"
#> [25] "TMe-2936" "TMe-3296" "TMe-3220" "TMe-3685" "TMe-2913" "TMe-3130"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2"    "TMe-2611" "TMe-3066" "TMe-85"   "TMe-1505" "TMe-2414"
#>  [7] "TMe-2935" "TMe-681"  "TMe-2123" "TMe-455"  "TMe-2242" "TMe-589" 
#> [13] "TMe-3187" "TMe-2705" "TMe-2843" "TMe-660"  "TMe-1833" "TMe-1864"
#> [19] "TMe-3188" "TMe-1831" "TMe-2352" "TMe-2077" "TMe-2021" "TMe-2568"
#> [25] "TMe-2851" "TMe-74"   "TMe-3093" "TMe-3605" "TMe-176"  "TMe-1422"
#> [31] "TMe-3338"
#> 
#> $III
#>  [1] "TMe-3207" "TMe-3028" "TMe-3166" "TMe-64"   "TMe-3324" "TMe-123" 
#>  [7] "TMe-1262" "TMe-3118" "TMe-161"  "TMe-1398" "TMe-206"  "TMe-1138"
#> [13] "TMe-3335" "TMe-2158" "TMe-2394" "TMe-1675" "TMe-3598" "TMe-3043"
#> 
#> $IV
#>  [1] "TMe-12"   "TMe-204"  "TMe-3107" "TMe-72"   "TMe-525"  "TMe-1364"
#>  [7] "TMe-2954" "TMe-1139" "TMe-1903" "TMe-2390" "TMe-3168" "TMe-317" 
#> [13] "TMe-3454" "TMe-2340" "TMe-552"  "TMe-396"  "TMe-1726" "TMe-2039"
#> [19] "TMe-1336" "TMe-3231" "TMe-2947" "TMe-885"  "TMe-840"  "TMe-190" 
#> [25] "TMe-802"  "TMe-2928" "TMe-1419" "TMe-3068" "TMe-3044" "TMe-2210"
#> [31] "TMe-2989" "TMe-3562" "TMe-3253" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-2361" "TMe-3332" "TMe-348"  "TMe-1127" "TMe-2400" "TMe-212" 
#>  [7] "TMe-1414" "TMe-436"  "TMe-479"  "TMe-1399" "TMe-1299" "TMe-359" 
#> [13] "TMe-826"  "TMe-1622" "TMe-1559" "TMe-168"  "TMe-1011" "TMe-685" 
#> [19] "TMe-417"  "TMe-1234" "TMe-424"  "TMe-667"  "TMe-1879" "TMe-456" 
#> [25] "TMe-629"  "TMe-982"  "TMe-1523" "TMe-1012" "TMe-797"  "TMe-669" 
#> [31] "TMe-1901" "TMe-2435" "TMe-2290" "TMe-2032" "TMe-1733" "TMe-1487"
#> [37] "TMe-2530" "TMe-2961" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-2195" "TMe-1832" "TMe-1110" "TMe-1512" "TMe-1403"
#>  [7] "TMe-809"  "TMe-1481" "TMe-2308" "TMe-2389" "TMe-1531" "TMe-3708"
#> [13] "TMe-1017" "TMe-281"  "TMe-791"  "TMe-625"  "TMe-1661"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_pam_out,
          use.names = FALSE)) +
  labs(title = "optim.medoid")


# Cluster-Based Sampling via K-means (Naes Method)
sel_naes_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "naes")
sel_naes_out
#> $I
#>  [1] "TMe-3026" "TMe-3282" "TMe-1532" "TMe-3087" "TMe-3076" "TMe-610" 
#>  [7] "TMe-306"  "TMe-2131" "TMe-1344" "TMe-2944" "TMe-1140" "TMe-486" 
#> [13] "TMe-3296" "TMe-3394" "TMe-3236" "TMe-2779" "TMe-3110" "TMe-3030"
#> [19] "TMe-3601" "TMe-1581" "TMe-2604" "TMe-606"  "TMe-1981" "TMe-28"  
#> [25] "TMe-3337" "TMe-765"  "TMe-3416" "TMe-3130" "TMe-1341" "TMe-1734"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-1619" "TMe-2304" "TMe-2951" "TMe-2428" "TMe-3805" "TMe-2765"
#>  [7] "TMe-3366" "TMe-40"   "TMe-3766" "TMe-2970" "TMe-3093" "TMe-1474"
#> [13] "TMe-3066" "TMe-539"  "TMe-1504" "TMe-1864" "TMe-1732" "TMe-2242"
#> [19] "TMe-660"  "TMe-2414" "TMe-1919" "TMe-404"  "TMe-377"  "TMe-589" 
#> [25] "TMe-2329" "TMe-1461" "TMe-2998" "TMe-2860" "TMe-2258" "TMe-1766"
#> [31] "TMe-455" 
#> 
#> $III
#>  [1] "TMe-13"   "TMe-1138" "TMe-1836" "TMe-1817" "TMe-1421" "TMe-3638"
#>  [7] "TMe-878"  "TMe-1863" "TMe-830"  "TMe-64"   "TMe-3094" "TMe-3750"
#> [13] "TMe-2926" "TMe-3046" "TMe-3088" "TMe-3005" "TMe-2811" "TMe-3804"
#> 
#> $IV
#>  [1] "TMe-241"  "TMe-2410" "TMe-2833" "TMe-675"  "TMe-552"  "TMe-3639"
#>  [7] "TMe-650"  "TMe-1364" "TMe-1336" "TMe-3269" "TMe-1402" "TMe-3568"
#> [13] "TMe-1485" "TMe-380"  "TMe-3168" "TMe-3535" "TMe-266"  "TMe-107" 
#> [19] "TMe-12"   "TMe-3000" "TMe-1765" "TMe-1297" "TMe-1139" "TMe-63"  
#> [25] "TMe-3253" "TMe-2971" "TMe-24"   "TMe-1377" "TMe-634"  "TMe-1867"
#> [31] "TMe-82"   "TMe-3585" "TMe-3242" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-645"  "TMe-2339" "TMe-3356" "TMe-769"  "TMe-181"  "TMe-363" 
#>  [7] "TMe-582"  "TMe-1188" "TMe-870"  "TMe-323"  "TMe-3363" "TMe-982" 
#> [13] "TMe-2518" "TMe-1924" "TMe-2769" "TMe-1345" "TMe-1293" "TMe-1859"
#> [19] "TMe-2530" "TMe-1227" "TMe-1788" "TMe-2688" "TMe-2904" "TMe-1979"
#> [25] "TMe-1388" "TMe-3329" "TMe-464"  "TMe-1572" "TMe-1099" "TMe-826" 
#> [31] "TMe-759"  "TMe-336"  "TMe-945"  "TMe-2439" "TMe-439"  "TMe-2761"
#> [37] "TMe-418"  "TMe-2933" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-361"  "TMe-2963" "TMe-781"  "TMe-1661" "TMe-1074" "TMe-281" 
#>  [7] "TMe-1233" "TMe-838"  "TMe-1539" "TMe-846"  "TMe-185"  "TMe-1985"
#> [13] "TMe-1076" "TMe-1858" "TMe-1481" "TMe-1403" "TMe-1392"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_naes_out,
          use.names = FALSE)) +
  labs(title = "naes")


# Cluster-Based Sampling via Hierarchical Clustering with Random Selection

## UPGMA
sel_hclust_random_out1 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random",
                  hclust.method = "average")
sel_hclust_random_out1
#> $I
#>  [1] "TMe-583"  "TMe-44"   "TMe-3002" "TMe-2785" "TMe-3424" "TMe-446" 
#>  [7] "TMe-410"  "TMe-882"  "TMe-500"  "TMe-566"  "TMe-3282" "TMe-306" 
#> [13] "TMe-815"  "TMe-2993" "TMe-3322" "TMe-1532" "TMe-1621" "TMe-3112"
#> [19] "TMe-2967" "TMe-3548" "TMe-3478" "TMe-3389" "TMe-3539" "TMe-3437"
#> [25] "TMe-3705" "TMe-2934" "TMe-3130" "TMe-3163" "TMe-3292" "TMe-3415"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-339"  "TMe-3284" "TMe-3209" "TMe-1184" "TMe-289"  "TMe-369" 
#>  [7] "TMe-1107" "TMe-2997" "TMe-1172" "TMe-1444" "TMe-3447" "TMe-902" 
#> [13] "TMe-1907" "TMe-171"  "TMe-2891" "TMe-1919" "TMe-2033" "TMe-2056"
#> [19] "TMe-3498" "TMe-3531" "TMe-47"   "TMe-3140" "TMe-3466" "TMe-3605"
#> [25] "TMe-3366" "TMe-2952" "TMe-3200" "TMe-3210" "TMe-3222" "TMe-3338"
#> [31] "TMe-3690"
#> 
#> $III
#>  [1] "TMe-3731" "TMe-3326" "TMe-3277" "TMe-830"  "TMe-1006" "TMe-878" 
#>  [7] "TMe-267"  "TMe-405"  "TMe-420"  "TMe-434"  "TMe-1787" "TMe-2205"
#> [13] "TMe-3274" "TMe-3230" "TMe-3804" "TMe-3299" "TMe-3551" "TMe-3596"
#> 
#> $IV
#>  [1] "TMe-12"   "TMe-3233" "TMe-3243" "TMe-1118" "TMe-3760" "TMe-967" 
#>  [7] "TMe-226"  "TMe-2839" "TMe-1597" "TMe-3097" "TMe-2954" "TMe-2922"
#> [13] "TMe-2043" "TMe-1828" "TMe-550"  "TMe-3542" "TMe-2247" "TMe-812" 
#> [19] "TMe-3144" "TMe-1342" "TMe-3054" "TMe-2960" "TMe-1434" "TMe-2809"
#> [25] "TMe-3256" "TMe-3417" "TMe-3273" "TMe-3581" "TMe-3602" "TMe-698" 
#> [31] "TMe-1020" "TMe-3730" "TMe-3442" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-1541" "TMe-532"  "TMe-2124" "TMe-1131" "TMe-294" 
#>  [7] "TMe-997"  "TMe-1523" "TMe-193"  "TMe-2518" "TMe-920"  "TMe-1099"
#> [13] "TMe-1199" "TMe-418"  "TMe-419"  "TMe-858"  "TMe-603"  "TMe-730" 
#> [19] "TMe-341"  "TMe-2355" "TMe-3408" "TMe-3411" "TMe-2531" "TMe-1501"
#> [25] "TMe-742"  "TMe-616"  "TMe-2121" "TMe-3034" "TMe-2435" "TMe-2589"
#> [31] "TMe-2790" "TMe-2853" "TMe-3356" "TMe-344"  "TMe-536"  "TMe-712" 
#> [37] "TMe-1042" "TMe-3736" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-836"  "TMe-514"  "TMe-2543" "TMe-1506" "TMe-1232" "TMe-188" 
#>  [7] "TMe-1646" "TMe-876"  "TMe-690"  "TMe-1539" "TMe-1383" "TMe-1124"
#> [13] "TMe-693"  "TMe-2983" "TMe-2791" "TMe-1769" "TMe-1174"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out1,
          use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "average")


## Single-linkage
sel_hclust_random_out2 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random",
                  hclust.method = "single")
sel_hclust_random_out2
#> $I
#>  [1] "TMe-3488" "TMe-41"   "TMe-500"  "TMe-566"  "TMe-815"  "TMe-867" 
#>  [7] "TMe-1466" "TMe-1564" "TMe-2027" "TMe-2773" "TMe-2785" "TMe-2967"
#> [13] "TMe-3151" "TMe-606"  "TMe-717"  "TMe-1218" "TMe-3229" "TMe-2934"
#> [19] "TMe-2984" "TMe-2985" "TMe-3112" "TMe-3163" "TMe-3292" "TMe-3389"
#> [25] "TMe-3415" "TMe-3437" "TMe-3471" "TMe-3514" "TMe-2940" "TMe-3478"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2935" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-404"  "TMe-768" 
#>  [7] "TMe-1385" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2056" "TMe-2204"
#> [13] "TMe-2851" "TMe-2891" "TMe-3009" "TMe-3140" "TMe-3766" "TMe-3466"
#> [19] "TMe-3605" "TMe-171"  "TMe-431"  "TMe-478"  "TMe-2352" "TMe-2568"
#> [25] "TMe-2952" "TMe-3200" "TMe-3210" "TMe-3222" "TMe-3338" "TMe-3690"
#> [31] "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-3274" "TMe-261"  "TMe-267"  "TMe-381" 
#>  [7] "TMe-785"  "TMe-1200" "TMe-3299" "TMe-3804" "TMe-35"   "TMe-174" 
#> [13] "TMe-635"  "TMe-2756" "TMe-3143" "TMe-3445" "TMe-3551" "TMe-3575"
#> 
#> $IV
#>  [1] "TMe-2552" "TMe-154"  "TMe-812"  "TMe-919"  "TMe-1376" "TMe-1434"
#>  [7] "TMe-2043" "TMe-2809" "TMe-2956" "TMe-2971" "TMe-3025" "TMe-3054"
#> [13] "TMe-3055" "TMe-3255" "TMe-3273" "TMe-3545" "TMe-180"  "TMe-372" 
#> [19] "TMe-663"  "TMe-698"  "TMe-737"  "TMe-1010" "TMe-1020" "TMe-1988"
#> [25] "TMe-2922" "TMe-3167" "TMe-3297" "TMe-3386" "TMe-3422" "TMe-3730"
#> [31] "TMe-3316" "TMe-3442" "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-920"  "TMe-378"  "TMe-418"  "TMe-861"  "TMe-1099"
#>  [7] "TMe-1366" "TMe-1458" "TMe-1482" "TMe-1963" "TMe-2050" "TMe-2121"
#> [13] "TMe-2531" "TMe-2589" "TMe-2688" "TMe-2790" "TMe-2853" "TMe-2863"
#> [19] "TMe-2973" "TMe-3185" "TMe-189"  "TMe-277"  "TMe-284"  "TMe-536" 
#> [25] "TMe-616"  "TMe-712"  "TMe-1305" "TMe-2003" "TMe-2139" "TMe-2355"
#> [31] "TMe-2905" "TMe-2953" "TMe-2959" "TMe-3034" "TMe-3311" "TMe-3356"
#> [37] "TMe-3408" "TMe-3736" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-281"  "TMe-465"  "TMe-1416" "TMe-1428" "TMe-1775" "TMe-1875"
#>  [7] "TMe-2064" "TMe-2534" "TMe-1053" "TMe-1124" "TMe-1174" "TMe-1232"
#> [13] "TMe-1769" "TMe-2510" "TMe-3095" "TMe-2983" "TMe-3116"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out2,
          use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "single")


## Complete-linkage
sel_hclust_random_out3 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random",
                  hclust.method = "complete")
sel_hclust_random_out3
#> $I
#>  [1] "TMe-1981" "TMe-264"  "TMe-3112" "TMe-2513" "TMe-1830" "TMe-1190"
#>  [7] "TMe-2031" "TMe-1564" "TMe-1096" "TMe-1466" "TMe-1964" "TMe-2909"
#> [13] "TMe-2934" "TMe-3392" "TMe-3145" "TMe-3264" "TMe-1533" "TMe-610" 
#> [19] "TMe-888"  "TMe-2069" "TMe-3488" "TMe-3726" "TMe-3115" "TMe-2916"
#> [25] "TMe-3425" "TMe-3296" "TMe-3151" "TMe-3437" "TMe-3685" "TMe-3462"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3150" "TMe-86"   "TMe-842"  "TMe-3009" "TMe-3557" "TMe-3239"
#>  [7] "TMe-478"  "TMe-3146" "TMe-3447" "TMe-2054" "TMe-477"  "TMe-589" 
#> [13] "TMe-3258" "TMe-1051" "TMe-3354" "TMe-2123" "TMe-2997" "TMe-1623"
#> [19] "TMe-1907" "TMe-2843" "TMe-176"  "TMe-1860" "TMe-1919" "TMe-2033"
#> [25] "TMe-2757" "TMe-2204" "TMe-74"   "TMe-3771" "TMe-47"   "TMe-3140"
#> [31] "TMe-2952"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-200"  "TMe-304"  "TMe-3230" "TMe-3133" "TMe-3029"
#>  [7] "TMe-1819" "TMe-10"   "TMe-1910" "TMe-1398" "TMe-635"  "TMe-434" 
#> [13] "TMe-2347" "TMe-3772" "TMe-3407" "TMe-3032" "TMe-3148" "TMe-23"  
#> 
#> $IV
#>  [1] "TMe-3044" "TMe-3422" "TMe-3206" "TMe-3780" "TMe-107"  "TMe-180" 
#>  [7] "TMe-1364" "TMe-2800" "TMe-824"  "TMe-266"  "TMe-276"  "TMe-225" 
#> [13] "TMe-1328" "TMe-415"  "TMe-1166" "TMe-386"  "TMe-396"  "TMe-3619"
#> [19] "TMe-2036" "TMe-3168" "TMe-3218" "TMe-3232" "TMe-184"  "TMe-575" 
#> [25] "TMe-78"   "TMe-3573" "TMe-3242" "TMe-2960" "TMe-2833" "TMe-2809"
#> [31] "TMe-3025" "TMe-3275" "TMe-3707" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-208"  "TMe-2907" "TMe-2821" "TMe-1305" "TMe-813"  "TMe-1268"
#>  [7] "TMe-1785" "TMe-382"  "TMe-1871" "TMe-2361" "TMe-334"  "TMe-1609"
#> [13] "TMe-436"  "TMe-1098" "TMe-1099" "TMe-1042" "TMe-2055" "TMe-1322"
#> [19] "TMe-2009" "TMe-924"  "TMe-1009" "TMe-1234" "TMe-623"  "TMe-782" 
#> [25] "TMe-832"  "TMe-997"  "TMe-861"  "TMe-2933" "TMe-1367" "TMe-612" 
#> [31] "TMe-1482" "TMe-1957" "TMe-1730" "TMe-1269" "TMe-3034" "TMe-2530"
#> [37] "TMe-2435" "TMe-284"  "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1756" "TMe-1985" "TMe-2534" "TMe-1477" "TMe-361"  "TMe-1403"
#>  [7] "TMe-2389" "TMe-1033" "TMe-580"  "TMe-1392" "TMe-983"  "TMe-2026"
#> [13] "TMe-2035" "TMe-1900" "TMe-1383" "TMe-1486" "TMe-893" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out3,
          use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "complete")


## Ward's D
sel_hclust_random_out4 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random",
                  hclust.method = "ward.D")
sel_hclust_random_out4
#> $I
#>  [1] "TMe-642"  "TMe-2773" "TMe-3325" "TMe-132"  "TMe-3424" "TMe-2513"
#>  [7] "TMe-3719" "TMe-3142" "TMe-1344" "TMe-3379" "TMe-566"  "TMe-1347"
#> [13] "TMe-3223" "TMe-1533" "TMe-3345" "TMe-3252" "TMe-990"  "TMe-1425"
#> [19] "TMe-865"  "TMe-1568" "TMe-1906" "TMe-3074" "TMe-2936" "TMe-3515"
#> [25] "TMe-727"  "TMe-3272" "TMe-3685" "TMe-2943" "TMe-3163" "TMe-2927"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2058" "TMe-39"   "TMe-3066" "TMe-3052" "TMe-2166" "TMe-2950"
#>  [7] "TMe-2257" "TMe-404"  "TMe-1795" "TMe-1969" "TMe-3203" "TMe-1616"
#> [13] "TMe-3187" "TMe-2957" "TMe-1279" "TMe-3200" "TMe-660"  "TMe-2200"
#> [19] "TMe-1005" "TMe-2258" "TMe-3365" "TMe-171"  "TMe-1766" "TMe-3140"
#> [25] "TMe-2568" "TMe-3498" "TMe-3338" "TMe-3455" "TMe-3612" "TMe-3540"
#> [31] "TMe-3286"
#> 
#> $III
#>  [1] "TMe-155"  "TMe-2356" "TMe-1817" "TMe-1868" "TMe-3302" "TMe-234" 
#>  [7] "TMe-3679" "TMe-853"  "TMe-3277" "TMe-267"  "TMe-420"  "TMe-3551"
#> [13] "TMe-2405" "TMe-3321" "TMe-3085" "TMe-3397" "TMe-2968" "TMe-3274"
#> 
#> $IV
#>  [1] "TMe-3206" "TMe-204"  "TMe-3291" "TMe-2971" "TMe-87"   "TMe-962" 
#>  [7] "TMe-226"  "TMe-575"  "TMe-714"  "TMe-180"  "TMe-3144" "TMe-320" 
#> [13] "TMe-1350" "TMe-2390" "TMe-766"  "TMe-3019" "TMe-403"  "TMe-184" 
#> [19] "TMe-656"  "TMe-526"  "TMe-675"  "TMe-3454" "TMe-1402" "TMe-1162"
#> [25] "TMe-1916" "TMe-186"  "TMe-3585" "TMe-3265" "TMe-2025" "TMe-2928"
#> [31] "TMe-3068" "TMe-3276" "TMe-2947" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-376"  "TMe-2905" "TMe-158"  "TMe-1534" "TMe-1295" "TMe-275" 
#>  [7] "TMe-2959" "TMe-2855" "TMe-1294" "TMe-2953" "TMe-1551" "TMe-208" 
#> [13] "TMe-197"  "TMe-798"  "TMe-162"  "TMe-373"  "TMe-1440" "TMe-448" 
#> [19] "TMe-1158" "TMe-800"  "TMe-764"  "TMe-1487" "TMe-1157" "TMe-429" 
#> [25] "TMe-655"  "TMe-668"  "TMe-926"  "TMe-755"  "TMe-612"  "TMe-2853"
#> [31] "TMe-326"  "TMe-2769" "TMe-2290" "TMe-927"  "TMe-2933" "TMe-1375"
#> [37] "TMe-2015" "TMe-194"  "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-652"  "TMe-1416" "TMe-269"  "TMe-1816" "TMe-1995" "TMe-1518"
#>  [7] "TMe-1531" "TMe-3177" "TMe-728"  "TMe-662"  "TMe-1847" "TMe-2818"
#> [13] "TMe-1875" "TMe-1261" "TMe-1178" "TMe-1744" "TMe-1554"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out4,
          use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "ward.D")


## WPGMA
sel_hclust_random_out5 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random",
                  hclust.method = "mcquitty")
sel_hclust_random_out5
#> $I
#>  [1] "TMe-1915" "TMe-41"   "TMe-117"  "TMe-132"  "TMe-566"  "TMe-1935"
#>  [7] "TMe-1218" "TMe-469"  "TMe-2027" "TMe-677"  "TMe-670"  "TMe-3252"
#> [13] "TMe-757"  "TMe-1266" "TMe-3330" "TMe-2372" "TMe-1532" "TMe-3112"
#> [19] "TMe-2967" "TMe-3111" "TMe-3433" "TMe-3219" "TMe-3272" "TMe-3518"
#> [25] "TMe-3685" "TMe-2943" "TMe-3475" "TMe-3163" "TMe-3292" "TMe-3415"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-47"   "TMe-2611" "TMe-3467" "TMe-3009" "TMe-1184" "TMe-431" 
#>  [7] "TMe-369"  "TMe-1241" "TMe-2048" "TMe-1461" "TMe-2329" "TMe-1005"
#> [13] "TMe-2903" "TMe-2257" "TMe-3286" "TMe-1616" "TMe-1831" "TMe-1860"
#> [19] "TMe-1907" "TMe-1919" "TMe-3805" "TMe-1271" "TMe-3498" "TMe-74"  
#> [25] "TMe-3140" "TMe-3466" "TMe-3605" "TMe-3200" "TMe-3530" "TMe-3338"
#> [31] "TMe-3800"
#> 
#> $III
#>  [1] "TMe-155"  "TMe-1398" "TMe-3700" "TMe-1147" "TMe-2897" "TMe-2203"
#>  [7] "TMe-878"  "TMe-425"  "TMe-635"  "TMe-405"  "TMe-2756" "TMe-785" 
#> [13] "TMe-3721" "TMe-3565" "TMe-116"  "TMe-3592" "TMe-3299" "TMe-3596"
#> 
#> $IV
#>  [1] "TMe-3206" "TMe-66"   "TMe-78"   "TMe-1325" "TMe-760"  "TMe-766" 
#>  [7] "TMe-215"  "TMe-2924" "TMe-684"  "TMe-1010" "TMe-3109" "TMe-30"  
#> [13] "TMe-186"  "TMe-656"  "TMe-1278" "TMe-1828" "TMe-1095" "TMe-87"  
#> [19] "TMe-550"  "TMe-1406" "TMe-3422" "TMe-1988" "TMe-3734" "TMe-3707"
#> [25] "TMe-622"  "TMe-1434" "TMe-3527" "TMe-2809" "TMe-3253" "TMe-3701"
#> [31] "TMe-27"   "TMe-698"  "TMe-3442" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-208"  "TMe-2060" "TMe-2003" "TMe-1877" "TMe-612"  "TMe-1098"
#>  [7] "TMe-294"  "TMe-447"  "TMe-334"  "TMe-362"  "TMe-527"  "TMe-700" 
#> [13] "TMe-287"  "TMe-1411" "TMe-1265" "TMe-723"  "TMe-2121" "TMe-2032"
#> [19] "TMe-162"  "TMe-729"  "TMe-826"  "TMe-2355" "TMe-2016" "TMe-2790"
#> [25] "TMe-1332" "TMe-2531" "TMe-194"  "TMe-1730" "TMe-2015" "TMe-2050"
#> [31] "TMe-3034" "TMe-2589" "TMe-2853" "TMe-3185" "TMe-344"  "TMe-877" 
#> [37] "TMe-2905" "TMe-2961" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1110" "TMe-514"  "TMe-2818" "TMe-781"  "TMe-3549" "TMe-465" 
#>  [7] "TMe-1079" "TMe-690"  "TMe-1809" "TMe-1416" "TMe-2510" "TMe-1875"
#> [13] "TMe-2983" "TMe-1053" "TMe-2791" "TMe-1319" "TMe-1174"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out5,
          use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "mcquitty")


## WPGMC
sel_hclust_random_out6 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random",
                  hclust.method = "median")
sel_hclust_random_out6
#> $I
#>  [1] "TMe-3323" "TMe-41"   "TMe-44"   "TMe-500"  "TMe-566"  "TMe-569" 
#>  [7] "TMe-765"  "TMe-1425" "TMe-1466" "TMe-1935" "TMe-2027" "TMe-2372"
#> [13] "TMe-2967" "TMe-3151" "TMe-3452" "TMe-3481" "TMe-717"  "TMe-1589"
#> [19] "TMe-2996" "TMe-3163" "TMe-3292" "TMe-3389" "TMe-3415" "TMe-3433"
#> [25] "TMe-3437" "TMe-3471" "TMe-3475" "TMe-3705" "TMe-2943" "TMe-3462"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2021" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-509"  "TMe-768" 
#>  [7] "TMe-1444" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2128" "TMe-2204"
#> [13] "TMe-3009" "TMe-3766" "TMe-3401" "TMe-3466" "TMe-3605" "TMe-47"  
#> [19] "TMe-1619" "TMe-2860" "TMe-2952" "TMe-3200" "TMe-3210" "TMe-3222"
#> [25] "TMe-3338" "TMe-3406" "TMe-3547" "TMe-3690" "TMe-3800" "TMe-3286"
#> [31] "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-1889" "TMe-104"  "TMe-141"  "TMe-425"  "TMe-785" 
#>  [7] "TMe-1200" "TMe-1993" "TMe-3299" "TMe-3804" "TMe-35"   "TMe-234" 
#> [13] "TMe-3596" "TMe-3143" "TMe-3274" "TMe-3397" "TMe-3445" "TMe-3565"
#> 
#> $IV
#>  [1] "TMe-3760" "TMe-154"  "TMe-427"  "TMe-550"  "TMe-812"  "TMe-1434"
#>  [7] "TMe-1525" "TMe-1806" "TMe-2043" "TMe-2755" "TMe-241"  "TMe-2956"
#> [13] "TMe-3576" "TMe-3025" "TMe-3054" "TMe-3055" "TMe-3089" "TMe-3255"
#> [19] "TMe-3256" "TMe-3273" "TMe-225"  "TMe-226"  "TMe-259"  "TMe-698" 
#> [25] "TMe-737"  "TMe-1020" "TMe-1511" "TMe-3189" "TMe-3701" "TMe-3730"
#> [31] "TMe-3053" "TMe-3442" "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-1188" "TMe-399"  "TMe-418"  "TMe-419"  "TMe-603" 
#>  [7] "TMe-755"  "TMe-1204" "TMe-1458" "TMe-2041" "TMe-2050" "TMe-2121"
#> [13] "TMe-2531" "TMe-2589" "TMe-2790" "TMe-2853" "TMe-2855" "TMe-2863"
#> [19] "TMe-2973" "TMe-3185" "TMe-481"  "TMe-536"  "TMe-616"  "TMe-712" 
#> [25] "TMe-723"  "TMe-816"  "TMe-1042" "TMe-1196" "TMe-1283" "TMe-2905"
#> [31] "TMe-2953" "TMe-3034" "TMe-3311" "TMe-3356" "TMe-3408" "TMe-3628"
#> [37] "TMe-3736" "TMe-3659" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1233" "TMe-269"  "TMe-531"  "TMe-725"  "TMe-1416" "TMe-693" 
#>  [7] "TMe-1875" "TMe-1995" "TMe-2026" "TMe-2983" "TMe-2534" "TMe-2791"
#> [13] "TMe-310"  "TMe-1053" "TMe-1174" "TMe-1769" "TMe-2510"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out6,
          use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "median")


## UPGMC
sel_hclust_random_out7 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random",
                  hclust.method = "centroid")
sel_hclust_random_out7
#> $I
#>  [1] "TMe-3465" "TMe-41"   "TMe-410"  "TMe-500"  "TMe-566"  "TMe-569" 
#>  [7] "TMe-815"  "TMe-867"  "TMe-1360" "TMe-1425" "TMe-1466" "TMe-2027"
#> [13] "TMe-2069" "TMe-2967" "TMe-717"  "TMe-2934" "TMe-2984" "TMe-2985"
#> [19] "TMe-3163" "TMe-3282" "TMe-3292" "TMe-3389" "TMe-3415" "TMe-3437"
#> [25] "TMe-3471" "TMe-3475" "TMe-3705" "TMe-2943" "TMe-3249" "TMe-3478"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-589"  "TMe-160"  "TMe-289"  "TMe-369"  "TMe-509"  "TMe-539" 
#>  [7] "TMe-902"  "TMe-1385" "TMe-1444" "TMe-1907" "TMe-1919" "TMe-2033"
#> [13] "TMe-2056" "TMe-2128" "TMe-2204" "TMe-2891" "TMe-3009" "TMe-3466"
#> [19] "TMe-3605" "TMe-171"  "TMe-478"  "TMe-1619" "TMe-2952" "TMe-3200"
#> [25] "TMe-3210" "TMe-3222" "TMe-3338" "TMe-3366" "TMe-3690" "TMe-3800"
#> [31] "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-2285" "TMe-261"  "TMe-267"  "TMe-381"  "TMe-425" 
#>  [7] "TMe-785"  "TMe-1868" "TMe-1897" "TMe-3299" "TMe-3804" "TMe-35"  
#> [13] "TMe-116"  "TMe-635"  "TMe-3596" "TMe-3274" "TMe-3445" "TMe-3565"
#> 
#> $IV
#>  [1] "TMe-1402" "TMe-154"  "TMe-427"  "TMe-650"  "TMe-812"  "TMe-1434"
#>  [7] "TMe-1525" "TMe-2043" "TMe-2971" "TMe-3025" "TMe-3054" "TMe-3055"
#> [13] "TMe-3089" "TMe-3255" "TMe-3256" "TMe-3273" "TMe-180"  "TMe-226" 
#> [19] "TMe-372"  "TMe-421"  "TMe-608"  "TMe-698"  "TMe-1020" "TMe-1129"
#> [25] "TMe-2567" "TMe-3189" "TMe-3297" "TMe-3417" "TMe-3701" "TMe-3707"
#> [31] "TMe-3730" "TMe-3442" "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-1862" "TMe-378"  "TMe-399"  "TMe-419"  "TMe-603"  "TMe-782" 
#>  [7] "TMe-861"  "TMe-1381" "TMe-1458" "TMe-1924" "TMe-2050" "TMe-2060"
#> [13] "TMe-2121" "TMe-2290" "TMe-2413" "TMe-2531" "TMe-2589" "TMe-2688"
#> [19] "TMe-2790" "TMe-2853" "TMe-2973" "TMe-3185" "TMe-189"  "TMe-536" 
#> [25] "TMe-616"  "TMe-712"  "TMe-723"  "TMe-1042" "TMe-2003" "TMe-2139"
#> [31] "TMe-2355" "TMe-2905" "TMe-2953" "TMe-3034" "TMe-3311" "TMe-3356"
#> [37] "TMe-3408" "TMe-3736" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1509" "TMe-465"  "TMe-725"  "TMe-1416" "TMe-1428" "TMe-1875"
#>  [7] "TMe-2035" "TMe-2064" "TMe-2534" "TMe-2791" "TMe-1124" "TMe-1174"
#> [13] "TMe-1232" "TMe-1769" "TMe-2510" "TMe-3387" "TMe-2983"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out7,
          use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "centroid")


## Ward's D2
sel_hclust_random_out8 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random",
                  hclust.method = "ward.D2")
sel_hclust_random_out8
#> $I
#>  [1] "TMe-841"  "TMe-44"   "TMe-1922" "TMe-2785" "TMe-2914" "TMe-2040"
#>  [7] "TMe-1564" "TMe-1823" "TMe-3394" "TMe-1272" "TMe-500"  "TMe-3170"
#> [13] "TMe-566"  "TMe-1091" "TMe-3127" "TMe-3314" "TMe-1532" "TMe-3479"
#> [19] "TMe-3266" "TMe-1425" "TMe-839"  "TMe-3074" "TMe-3319" "TMe-3110"
#> [25] "TMe-3462" "TMe-3577" "TMe-3437" "TMe-3705" "TMe-2943" "TMe-3461"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3150" "TMe-2611" "TMe-842"  "TMe-1250" "TMe-2414" "TMe-339" 
#>  [7] "TMe-681"  "TMe-1795" "TMe-1969" "TMe-890"  "TMe-2329" "TMe-706" 
#> [13] "TMe-2951" "TMe-3531" "TMe-601"  "TMe-3200" "TMe-2866" "TMe-1005"
#> [19] "TMe-1241" "TMe-1860" "TMe-3447" "TMe-171"  "TMe-3498" "TMe-2568"
#> [25] "TMe-2033" "TMe-2860" "TMe-2915" "TMe-3671" "TMe-3020" "TMe-3140"
#> [31] "TMe-3338"
#> 
#> $III
#>  [1] "TMe-10"   "TMe-3029" "TMe-2977" "TMe-304"  "TMe-3572" "TMe-3772"
#>  [7] "TMe-148"  "TMe-1443" "TMe-206"  "TMe-1993" "TMe-3407" "TMe-2158"
#> [13] "TMe-1230" "TMe-1939" "TMe-3118" "TMe-35"   "TMe-3176" "TMe-3043"
#> 
#> $IV
#>  [1] "TMe-25"   "TMe-3606" "TMe-1456" "TMe-2839" "TMe-3055" "TMe-3273"
#>  [7] "TMe-270"  "TMe-1525" "TMe-460"  "TMe-2552" "TMe-766"  "TMe-735" 
#> [13] "TMe-3760" "TMe-3541" "TMe-2390" "TMe-3212" "TMe-396"  "TMe-1867"
#> [19] "TMe-428"  "TMe-3168" "TMe-3206" "TMe-641"  "TMe-2247" "TMe-3573"
#> [25] "TMe-3576" "TMe-43"   "TMe-3225" "TMe-3256" "TMe-3537" "TMe-421" 
#> [31] "TMe-2020" "TMe-3290" "TMe-3275" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-2643" "TMe-348"  "TMe-158"  "TMe-2790" "TMe-1500" "TMe-1042"
#>  [7] "TMe-2750" "TMe-479"  "TMe-543"  "TMe-582"  "TMe-336"  "TMe-224" 
#> [13] "TMe-390"  "TMe-414"  "TMe-1098" "TMe-764"  "TMe-378"  "TMe-863" 
#> [19] "TMe-402"  "TMe-800"  "TMe-723"  "TMe-1458" "TMe-2413" "TMe-448" 
#> [25] "TMe-623"  "TMe-592"  "TMe-730"  "TMe-341"  "TMe-2853" "TMe-2290"
#> [31] "TMe-647"  "TMe-2961" "TMe-1367" "TMe-1730" "TMe-1159" "TMe-2119"
#> [37] "TMe-3628" "TMe-2590" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-854"  "TMe-3387" "TMe-1017" "TMe-505"  "TMe-1995" "TMe-620" 
#>  [7] "TMe-1416" "TMe-1992" "TMe-631"  "TMe-781"  "TMe-752"  "TMe-1074"
#> [13] "TMe-1858" "TMe-1261" "TMe-693"  "TMe-1744" "TMe-1539"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out8,
          use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "ward.D2")


# Cluster-Based Sampling via Hierarchical Clustering with Medoid Selection

## UPGMA
sel_hclust_medoid_out1 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid",
                  hclust.method = "average")
sel_hclust_medoid_out1
#> $I
#>  [1] "TMe-1981" "TMe-2810" "TMe-1696" "TMe-300"  "TMe-3424" "TMe-1823"
#>  [7] "TMe-1564" "TMe-882"  "TMe-500"  "TMe-3441" "TMe-2604" "TMe-306" 
#> [13] "TMe-815"  "TMe-3330" "TMe-610"  "TMe-2253" "TMe-3490" "TMe-3026"
#> [19] "TMe-2967" "TMe-2964" "TMe-3268" "TMe-3151" "TMe-3419" "TMe-2985"
#> [25] "TMe-3685" "TMe-2927" "TMe-2996" "TMe-3163" "TMe-3292" "TMe-3415"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3617" "TMe-3066" "TMe-1250" "TMe-2414" "TMe-289"  "TMe-369" 
#>  [7] "TMe-2123" "TMe-1969" "TMe-3203" "TMe-2329" "TMe-2951" "TMe-902" 
#> [13] "TMe-674"  "TMe-176"  "TMe-2021" "TMe-1919" "TMe-2033" "TMe-2056"
#> [19] "TMe-3498" "TMe-3642" "TMe-3101" "TMe-3766" "TMe-3466" "TMe-3605"
#> [25] "TMe-1668" "TMe-2952" "TMe-3200" "TMe-3210" "TMe-3222" "TMe-3338"
#> [31] "TMe-3690"
#> 
#> $III
#>  [1] "TMe-93"   "TMe-3326" "TMe-3133" "TMe-123"  "TMe-187"  "TMe-1138"
#>  [7] "TMe-267"  "TMe-405"  "TMe-3638" "TMe-434"  "TMe-1787" "TMe-1288"
#> [13] "TMe-3663" "TMe-3592" "TMe-3721" "TMe-3299" "TMe-3598" "TMe-2926"
#> 
#> $IV
#>  [1] "TMe-3238" "TMe-1948" "TMe-3780" "TMe-1470" "TMe-81"   "TMe-108" 
#>  [7] "TMe-802"  "TMe-2928" "TMe-684"  "TMe-314"  "TMe-2954" "TMe-2390"
#> [13] "TMe-396"  "TMe-1336" "TMe-550"  "TMe-3542" "TMe-1139" "TMe-812" 
#> [19] "TMe-897"  "TMe-2567" "TMe-1350" "TMe-3585" "TMe-1434" "TMe-2809"
#> [25] "TMe-2971" "TMe-3591" "TMe-3273" "TMe-3581" "TMe-3602" "TMe-698" 
#> [31] "TMe-1020" "TMe-3730" "TMe-3442" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-99"   "TMe-1127" "TMe-2009" "TMe-645"  "TMe-294" 
#>  [7] "TMe-332"  "TMe-582"  "TMe-2643" "TMe-826"  "TMe-1234" "TMe-1099"
#> [13] "TMe-800"  "TMe-418"  "TMe-1487" "TMe-2032" "TMe-603"  "TMe-730" 
#> [19] "TMe-341"  "TMe-861"  "TMe-813"  "TMe-2933" "TMe-1307" "TMe-1957"
#> [25] "TMe-2753" "TMe-2050" "TMe-2121" "TMe-2290" "TMe-2435" "TMe-2589"
#> [31] "TMe-2790" "TMe-2853" "TMe-3185" "TMe-284"  "TMe-536"  "TMe-712" 
#> [37] "TMe-1042" "TMe-2905" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-2308" "TMe-791"  "TMe-1403" "TMe-465"  "TMe-281" 
#>  [7] "TMe-2551" "TMe-1216" "TMe-2963" "TMe-2389" "TMe-1383" "TMe-1416"
#> [13] "TMe-1445" "TMe-2035" "TMe-2791" "TMe-1053" "TMe-1174"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out1,
          use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "average")


## Single-linkage
sel_hclust_medoid_out2 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid",
                  hclust.method = "single")
sel_hclust_medoid_out2
#> $I
#>  [1] "TMe-1696" "TMe-41"   "TMe-500"  "TMe-566"  "TMe-815"  "TMe-867" 
#>  [7] "TMe-1466" "TMe-1564" "TMe-2027" "TMe-2773" "TMe-2785" "TMe-2967"
#> [13] "TMe-3151" "TMe-606"  "TMe-717"  "TMe-1218" "TMe-2462" "TMe-2934"
#> [19] "TMe-2984" "TMe-2985" "TMe-3112" "TMe-3163" "TMe-3292" "TMe-3389"
#> [25] "TMe-3415" "TMe-3437" "TMe-3471" "TMe-3514" "TMe-2940" "TMe-3478"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-404"  "TMe-768" 
#>  [7] "TMe-1385" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2056" "TMe-2204"
#> [13] "TMe-2851" "TMe-2891" "TMe-3009" "TMe-3140" "TMe-3141" "TMe-3466"
#> [19] "TMe-3605" "TMe-176"  "TMe-431"  "TMe-478"  "TMe-2352" "TMe-2568"
#> [25] "TMe-2952" "TMe-3200" "TMe-3210" "TMe-3222" "TMe-3338" "TMe-3690"
#> [31] "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-2154" "TMe-261"  "TMe-267"  "TMe-381" 
#>  [7] "TMe-785"  "TMe-1200" "TMe-3299" "TMe-3804" "TMe-35"   "TMe-174" 
#> [13] "TMe-635"  "TMe-2756" "TMe-3143" "TMe-3445" "TMe-3551" "TMe-3575"
#> 
#> $IV
#>  [1] "TMe-2340" "TMe-154"  "TMe-812"  "TMe-919"  "TMe-1376" "TMe-1434"
#>  [7] "TMe-2043" "TMe-2809" "TMe-2956" "TMe-2971" "TMe-3025" "TMe-3054"
#> [13] "TMe-3055" "TMe-3255" "TMe-3273" "TMe-3545" "TMe-180"  "TMe-372" 
#> [19] "TMe-663"  "TMe-698"  "TMe-737"  "TMe-1010" "TMe-1020" "TMe-1988"
#> [25] "TMe-2922" "TMe-3167" "TMe-3297" "TMe-3386" "TMe-3422" "TMe-3730"
#> [31] "TMe-3316" "TMe-3442" "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-332"  "TMe-378"  "TMe-418"  "TMe-861"  "TMe-1099"
#>  [7] "TMe-1366" "TMe-1458" "TMe-1482" "TMe-1963" "TMe-2050" "TMe-2121"
#> [13] "TMe-2531" "TMe-2589" "TMe-2688" "TMe-2790" "TMe-2853" "TMe-2863"
#> [19] "TMe-2973" "TMe-3185" "TMe-189"  "TMe-277"  "TMe-284"  "TMe-536" 
#> [25] "TMe-616"  "TMe-712"  "TMe-1305" "TMe-2003" "TMe-2139" "TMe-2355"
#> [31] "TMe-2905" "TMe-2953" "TMe-2959" "TMe-3034" "TMe-3311" "TMe-3356"
#> [37] "TMe-3408" "TMe-3736" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-465"  "TMe-1416" "TMe-1428" "TMe-1775" "TMe-1875"
#>  [7] "TMe-2064" "TMe-2534" "TMe-1053" "TMe-1124" "TMe-1174" "TMe-1232"
#> [13] "TMe-1769" "TMe-2510" "TMe-3095" "TMe-2983" "TMe-3116"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out2,
          use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "single")


## Complete-linkage
sel_hclust_medoid_out3 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid",
                  hclust.method = "complete")
sel_hclust_medoid_out3
#> $I
#>  [1] "TMe-1981" "TMe-1621" "TMe-3490" "TMe-1696" "TMe-3480" "TMe-3424"
#>  [7] "TMe-1930" "TMe-1564" "TMe-882"  "TMe-306"  "TMe-3322" "TMe-3496"
#> [13] "TMe-3051" "TMe-2604" "TMe-28"   "TMe-3330" "TMe-1086" "TMe-610" 
#> [19] "TMe-1272" "TMe-1532" "TMe-3518" "TMe-2936" "TMe-3065" "TMe-2913"
#> [25] "TMe-3425" "TMe-3235" "TMe-3151" "TMe-3389" "TMe-3667" "TMe-3104"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2"    "TMe-3533" "TMe-3571" "TMe-3009" "TMe-2414" "TMe-377" 
#>  [7] "TMe-2715" "TMe-2765" "TMe-433"  "TMe-1833" "TMe-890"  "TMe-589" 
#> [13] "TMe-2951" "TMe-601"  "TMe-3564" "TMe-2268" "TMe-3365" "TMe-2048"
#> [19] "TMe-674"  "TMe-1250" "TMe-176"  "TMe-2021" "TMe-1668" "TMe-2033"
#> [25] "TMe-2995" "TMe-3498" "TMe-3642" "TMe-3093" "TMe-3101" "TMe-3766"
#> [31] "TMe-2952"
#> 
#> $III
#>  [1] "TMe-155"  "TMe-1443" "TMe-3028" "TMe-93"   "TMe-617"  "TMe-1421"
#>  [7] "TMe-1262" "TMe-161"  "TMe-1138" "TMe-3321" "TMe-3362" "TMe-434" 
#> [13] "TMe-2394" "TMe-3094" "TMe-3324" "TMe-3750" "TMe-3043" "TMe-23"  
#> 
#> $IV
#>  [1] "TMe-3073" "TMe-3045" "TMe-2946" "TMe-3000" "TMe-3265" "TMe-634" 
#>  [7] "TMe-1364" "TMe-2410" "TMe-1027" "TMe-1903" "TMe-65"   "TMe-3761"
#> [13] "TMe-1419" "TMe-387"  "TMe-3541" "TMe-386"  "TMe-396"  "TMe-1726"
#> [19] "TMe-1336" "TMe-802"  "TMe-1470" "TMe-1231" "TMe-550"  "TMe-3108"
#> [25] "TMe-1144" "TMe-3019" "TMe-3068" "TMe-3585" "TMe-2059" "TMe-2809"
#> [31] "TMe-3025" "TMe-3253" "TMe-3255" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-1877" "TMe-98"   "TMe-1127" "TMe-1305" "TMe-212"  "TMe-645" 
#>  [7] "TMe-294"  "TMe-629"  "TMe-307"  "TMe-2441" "TMe-1399" "TMe-2213"
#> [13] "TMe-1358" "TMe-982"  "TMe-1455" "TMe-1011" "TMe-402"  "TMe-424" 
#> [19] "TMe-2009" "TMe-870"  "TMe-1299" "TMe-1234" "TMe-623"  "TMe-730" 
#> [25] "TMe-1227" "TMe-997"  "TMe-2355" "TMe-2933" "TMe-1901" "TMe-137" 
#> [31] "TMe-1963" "TMe-1500" "TMe-685"  "TMe-1733" "TMe-3034" "TMe-2530"
#> [37] "TMe-2435" "TMe-284"  "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-1985" "TMe-1164" "TMe-1721" "TMe-548"  "TMe-620" 
#>  [7] "TMe-1633" "TMe-838"  "TMe-1816" "TMe-2402" "TMe-281"  "TMe-625" 
#> [13] "TMe-2963" "TMe-1178" "TMe-3095" "TMe-2818" "TMe-1646"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out3,
          use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "complete")


## Ward's D
sel_hclust_medoid_out4 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid",
                  hclust.method = "ward.D")
sel_hclust_medoid_out4
#> $I
#>  [1] "TMe-2779" "TMe-2810" "TMe-1717" "TMe-300"  "TMe-3424" "TMe-1140"
#>  [7] "TMe-3465" "TMe-489"  "TMe-3202" "TMe-3548" "TMe-3441" "TMe-1981"
#> [13] "TMe-2604" "TMe-2253" "TMe-3132" "TMe-3252" "TMe-1823" "TMe-610" 
#> [19] "TMe-865"  "TMe-896"  "TMe-486"  "TMe-3477" "TMe-2936" "TMe-3268"
#> [25] "TMe-3220" "TMe-2985" "TMe-3685" "TMe-3694" "TMe-2964" "TMe-3425"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-40"   "TMe-2611" "TMe-3571" "TMe-85"   "TMe-2414" "TMe-377" 
#>  [7] "TMe-681"  "TMe-478"  "TMe-2765" "TMe-1474" "TMe-455"  "TMe-706" 
#> [13] "TMe-2951" "TMe-196"  "TMe-1484" "TMe-2843" "TMe-2824" "TMe-2200"
#> [19] "TMe-674"  "TMe-2021" "TMe-1698" "TMe-176"  "TMe-1833" "TMe-2352"
#> [25] "TMe-2056" "TMe-3308" "TMe-74"   "TMe-3455" "TMe-3093" "TMe-3540"
#> [31] "TMe-3188"
#> 
#> $III
#>  [1] "TMe-2161" "TMe-1288" "TMe-3592" "TMe-1836" "TMe-3302" "TMe-187" 
#>  [7] "TMe-206"  "TMe-1138" "TMe-3335" "TMe-3362" "TMe-2158" "TMe-3598"
#> [13] "TMe-2394" "TMe-3321" "TMe-3118" "TMe-3324" "TMe-3207" "TMe-3750"
#> 
#> $IV
#>  [1] "TMe-12"   "TMe-204"  "TMe-62"   "TMe-3232" "TMe-107"  "TMe-634" 
#>  [7] "TMe-1364" "TMe-3107" "TMe-1027" "TMe-415"  "TMe-3191" "TMe-2989"
#> [13] "TMe-3541" "TMe-386"  "TMe-387"  "TMe-396"  "TMe-1651" "TMe-3168"
#> [19] "TMe-885"  "TMe-2039" "TMe-2332" "TMe-891"  "TMe-1402" "TMe-1419"
#> [25] "TMe-2399" "TMe-1485" "TMe-3585" "TMe-1336" "TMe-1903" "TMe-2928"
#> [31] "TMe-3068" "TMe-3253" "TMe-63"   "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-804"  "TMe-307"  "TMe-1127" "TMe-1694" "TMe-212"  "TMe-868" 
#>  [7] "TMe-813"  "TMe-510"  "TMe-1160" "TMe-1399" "TMe-1879" "TMe-2643"
#> [13] "TMe-1358" "TMe-826"  "TMe-1234" "TMe-629"  "TMe-1730" "TMe-1100"
#> [19] "TMe-1011" "TMe-800"  "TMe-217"  "TMe-1487" "TMe-2032" "TMe-429" 
#> [25] "TMe-456"  "TMe-945"  "TMe-2845" "TMe-730"  "TMe-612"  "TMe-1307"
#> [31] "TMe-1411" "TMe-1760" "TMe-2290" "TMe-1188" "TMe-2933" "TMe-1901"
#> [37] "TMe-1733" "TMe-2435" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-2308" "TMe-1017" "TMe-1217" "TMe-907"  "TMe-1403"
#>  [7] "TMe-838"  "TMe-548"  "TMe-705"  "TMe-281"  "TMe-1164" "TMe-1835"
#> [13] "TMe-1079" "TMe-2963" "TMe-625"  "TMe-1661" "TMe-2389"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out4,
          use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "ward.D")


## WPGMA
sel_hclust_medoid_out5 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid",
                  hclust.method = "mcquitty")
sel_hclust_medoid_out5
#> $I
#>  [1] "TMe-1981" "TMe-2810" "TMe-117"  "TMe-300"  "TMe-3424" "TMe-1226"
#>  [7] "TMe-1564" "TMe-3202" "TMe-500"  "TMe-2964" "TMe-910"  "TMe-2604"
#> [13] "TMe-3425" "TMe-2066" "TMe-3330" "TMe-610"  "TMe-3490" "TMe-3493"
#> [19] "TMe-2967" "TMe-3268" "TMe-3151" "TMe-3419" "TMe-2985" "TMe-2917"
#> [25] "TMe-3685" "TMe-3130" "TMe-2927" "TMe-3163" "TMe-3292" "TMe-3415"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3540" "TMe-3533" "TMe-3571" "TMe-160"  "TMe-2414" "TMe-431" 
#>  [7] "TMe-369"  "TMe-1107" "TMe-1864" "TMe-1833" "TMe-2329" "TMe-589" 
#> [13] "TMe-2824" "TMe-2797" "TMe-1827" "TMe-1422" "TMe-1831" "TMe-2021"
#> [19] "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2568" "TMe-3498" "TMe-3642"
#> [25] "TMe-3766" "TMe-3466" "TMe-3605" "TMe-3200" "TMe-3284" "TMe-3338"
#> [31] "TMe-3800"
#> 
#> $III
#>  [1] "TMe-2161" "TMe-1443" "TMe-3069" "TMe-93"   "TMe-1819" "TMe-2203"
#>  [7] "TMe-1138" "TMe-3638" "TMe-267"  "TMe-405"  "TMe-3028" "TMe-785" 
#> [13] "TMe-3721" "TMe-3094" "TMe-3750" "TMe-3592" "TMe-3299" "TMe-2926"
#> 
#> $IV
#>  [1] "TMe-3454" "TMe-3000" "TMe-3068" "TMe-1652" "TMe-3265" "TMe-1094"
#>  [7] "TMe-1364" "TMe-2928" "TMe-1155" "TMe-415"  "TMe-1148" "TMe-2954"
#> [13] "TMe-1211" "TMe-2807" "TMe-2390" "TMe-3639" "TMe-3168" "TMe-540" 
#> [19] "TMe-550"  "TMe-812"  "TMe-919"  "TMe-1377" "TMe-2980" "TMe-3585"
#> [25] "TMe-2947" "TMe-1434" "TMe-153"  "TMe-2809" "TMe-3591" "TMe-3255"
#> [31] "TMe-3581" "TMe-698"  "TMe-3442" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-3332" "TMe-1871" "TMe-1127" "TMe-1877" "TMe-1722" "TMe-332" 
#>  [7] "TMe-197"  "TMe-456"  "TMe-1399" "TMe-362"  "TMe-1011" "TMe-1099"
#> [13] "TMe-629"  "TMe-245"  "TMe-800"  "TMe-418"  "TMe-1455" "TMe-2032"
#> [19] "TMe-1234" "TMe-730"  "TMe-2530" "TMe-861"  "TMe-1367" "TMe-972" 
#> [25] "TMe-1357" "TMe-1307" "TMe-2435" "TMe-2753" "TMe-1269" "TMe-2050"
#> [31] "TMe-2290" "TMe-2589" "TMe-2853" "TMe-3185" "TMe-284"  "TMe-341" 
#> [37] "TMe-2905" "TMe-2961" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-2308" "TMe-1074" "TMe-907"  "TMe-1403" "TMe-465" 
#>  [7] "TMe-625"  "TMe-2963" "TMe-620"  "TMe-1416" "TMe-2389" "TMe-2551"
#> [13] "TMe-2035" "TMe-1608" "TMe-2791" "TMe-188"  "TMe-1174"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out5,
          use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "mcquitty")


## WPGMC
sel_hclust_medoid_out6 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid",
                  hclust.method = "median")
sel_hclust_medoid_out6
#> $I
#>  [1] "TMe-1696" "TMe-41"   "TMe-44"   "TMe-500"  "TMe-566"  "TMe-569" 
#>  [7] "TMe-1190" "TMe-1425" "TMe-1466" "TMe-1935" "TMe-2027" "TMe-2372"
#> [13] "TMe-2967" "TMe-3151" "TMe-3452" "TMe-3481" "TMe-717"  "TMe-1589"
#> [19] "TMe-2996" "TMe-3163" "TMe-3292" "TMe-3389" "TMe-3415" "TMe-3433"
#> [25] "TMe-3437" "TMe-3471" "TMe-3475" "TMe-3705" "TMe-2943" "TMe-3462"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-509"  "TMe-768" 
#>  [7] "TMe-1444" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2128" "TMe-2204"
#> [13] "TMe-3009" "TMe-3766" "TMe-3401" "TMe-3466" "TMe-3605" "TMe-47"  
#> [19] "TMe-1619" "TMe-2860" "TMe-2952" "TMe-3200" "TMe-3210" "TMe-3222"
#> [25] "TMe-3338" "TMe-3406" "TMe-3547" "TMe-3690" "TMe-3800" "TMe-3286"
#> [31] "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-2154" "TMe-104"  "TMe-141"  "TMe-425"  "TMe-785" 
#>  [7] "TMe-1200" "TMe-1993" "TMe-3299" "TMe-3804" "TMe-35"   "TMe-234" 
#> [13] "TMe-2926" "TMe-3143" "TMe-3274" "TMe-3397" "TMe-3445" "TMe-3565"
#> 
#> $IV
#>  [1] "TMe-1726" "TMe-154"  "TMe-427"  "TMe-550"  "TMe-812"  "TMe-1434"
#>  [7] "TMe-1525" "TMe-1806" "TMe-2043" "TMe-2755" "TMe-2809" "TMe-2956"
#> [13] "TMe-3576" "TMe-3025" "TMe-3054" "TMe-3055" "TMe-3089" "TMe-3255"
#> [19] "TMe-3256" "TMe-3273" "TMe-225"  "TMe-226"  "TMe-259"  "TMe-698" 
#> [25] "TMe-737"  "TMe-1020" "TMe-1511" "TMe-3189" "TMe-3701" "TMe-3730"
#> [31] "TMe-3053" "TMe-3442" "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-332"  "TMe-399"  "TMe-418"  "TMe-419"  "TMe-603" 
#>  [7] "TMe-755"  "TMe-828"  "TMe-1458" "TMe-2041" "TMe-2050" "TMe-2121"
#> [13] "TMe-2531" "TMe-2589" "TMe-2790" "TMe-2853" "TMe-2855" "TMe-2863"
#> [19] "TMe-2973" "TMe-3185" "TMe-481"  "TMe-536"  "TMe-616"  "TMe-712" 
#> [25] "TMe-723"  "TMe-816"  "TMe-1042" "TMe-1196" "TMe-1283" "TMe-2905"
#> [31] "TMe-2953" "TMe-3034" "TMe-3311" "TMe-3356" "TMe-3408" "TMe-3628"
#> [37] "TMe-3736" "TMe-3659" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-269"  "TMe-531"  "TMe-725"  "TMe-1416" "TMe-1445"
#>  [7] "TMe-1875" "TMe-1995" "TMe-2026" "TMe-2035" "TMe-2534" "TMe-2791"
#> [13] "TMe-310"  "TMe-1053" "TMe-1174" "TMe-1769" "TMe-2510"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out6,
          use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "median")


## UPGMC
sel_hclust_medoid_out7 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid",
                  hclust.method = "centroid")
sel_hclust_medoid_out7
#> $I
#>  [1] "TMe-1696" "TMe-41"   "TMe-410"  "TMe-500"  "TMe-566"  "TMe-569" 
#>  [7] "TMe-815"  "TMe-867"  "TMe-1360" "TMe-1425" "TMe-1466" "TMe-2027"
#> [13] "TMe-2069" "TMe-2967" "TMe-717"  "TMe-2934" "TMe-2984" "TMe-2985"
#> [19] "TMe-3163" "TMe-3282" "TMe-3292" "TMe-3389" "TMe-3415" "TMe-3437"
#> [25] "TMe-3471" "TMe-3475" "TMe-3705" "TMe-2943" "TMe-3249" "TMe-3478"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-509"  "TMe-539" 
#>  [7] "TMe-902"  "TMe-1385" "TMe-1444" "TMe-1907" "TMe-1919" "TMe-2033"
#> [13] "TMe-2056" "TMe-2128" "TMe-2204" "TMe-2891" "TMe-3009" "TMe-3466"
#> [19] "TMe-3605" "TMe-171"  "TMe-478"  "TMe-1619" "TMe-2952" "TMe-3200"
#> [25] "TMe-3210" "TMe-3222" "TMe-3338" "TMe-3366" "TMe-3690" "TMe-3800"
#> [31] "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-2154" "TMe-261"  "TMe-267"  "TMe-381"  "TMe-425" 
#>  [7] "TMe-785"  "TMe-1868" "TMe-1897" "TMe-3299" "TMe-3804" "TMe-35"  
#> [13] "TMe-116"  "TMe-635"  "TMe-2926" "TMe-3274" "TMe-3445" "TMe-3565"
#> 
#> $IV
#>  [1] "TMe-2340" "TMe-154"  "TMe-427"  "TMe-650"  "TMe-812"  "TMe-1434"
#>  [7] "TMe-1525" "TMe-2043" "TMe-2971" "TMe-3025" "TMe-3054" "TMe-3055"
#> [13] "TMe-3089" "TMe-3255" "TMe-3256" "TMe-3273" "TMe-180"  "TMe-226" 
#> [19] "TMe-372"  "TMe-421"  "TMe-608"  "TMe-698"  "TMe-1020" "TMe-1129"
#> [25] "TMe-2567" "TMe-3189" "TMe-3297" "TMe-3417" "TMe-3701" "TMe-3707"
#> [31] "TMe-3730" "TMe-3442" "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-378"  "TMe-399"  "TMe-419"  "TMe-603"  "TMe-782" 
#>  [7] "TMe-861"  "TMe-1381" "TMe-1458" "TMe-1924" "TMe-2050" "TMe-2060"
#> [13] "TMe-2121" "TMe-2290" "TMe-2413" "TMe-2531" "TMe-2589" "TMe-2688"
#> [19] "TMe-2790" "TMe-2853" "TMe-2973" "TMe-3185" "TMe-189"  "TMe-536" 
#> [25] "TMe-616"  "TMe-712"  "TMe-723"  "TMe-1042" "TMe-2003" "TMe-2139"
#> [31] "TMe-2355" "TMe-2905" "TMe-2953" "TMe-3034" "TMe-3311" "TMe-3356"
#> [37] "TMe-3408" "TMe-3736" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-465"  "TMe-725"  "TMe-1416" "TMe-1428" "TMe-1875"
#>  [7] "TMe-2035" "TMe-2064" "TMe-2534" "TMe-2791" "TMe-1124" "TMe-1174"
#> [13] "TMe-1232" "TMe-1769" "TMe-2510" "TMe-3387" "TMe-2983"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out7,
          use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "centroid")


## Ward's D2
sel_hclust_medoid_out8 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts,
                  dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid",
                  hclust.method = "ward.D2")
sel_hclust_medoid_out8
#> $I
#>  [1] "TMe-642"  "TMe-2810" "TMe-1717" "TMe-300"  "TMe-3424" "TMe-486" 
#>  [7] "TMe-3465" "TMe-1823" "TMe-3202" "TMe-882"  "TMe-306"  "TMe-2964"
#> [13] "TMe-3441" "TMe-1448" "TMe-1981" "TMe-756"  "TMe-2253" "TMe-3132"
#> [19] "TMe-3252" "TMe-610"  "TMe-1086" "TMe-3477" "TMe-2604" "TMe-2936"
#> [25] "TMe-3419" "TMe-3577" "TMe-2985" "TMe-3685" "TMe-3130" "TMe-3425"
#> [31] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-40"   "TMe-2611" "TMe-3571" "TMe-85"   "TMe-2414" "TMe-377" 
#>  [7] "TMe-2715" "TMe-2765" "TMe-1474" "TMe-3203" "TMe-2329" "TMe-706" 
#> [13] "TMe-2951" "TMe-196"  "TMe-601"  "TMe-2843" "TMe-2824" "TMe-674" 
#> [19] "TMe-2123" "TMe-2021" "TMe-1698" "TMe-176"  "TMe-1833" "TMe-2568"
#> [25] "TMe-2033" "TMe-74"   "TMe-3455" "TMe-3093" "TMe-3540" "TMe-3766"
#> [31] "TMe-2952"
#> 
#> $III
#>  [1] "TMe-3362" "TMe-1288" "TMe-3592" "TMe-1836" "TMe-3302" "TMe-3324"
#>  [7] "TMe-187"  "TMe-3028" "TMe-206"  "TMe-1138" "TMe-3005" "TMe-2158"
#> [13] "TMe-2394" "TMe-3321" "TMe-3118" "TMe-3207" "TMe-3750" "TMe-3043"
#> 
#> $IV
#>  [1] "TMe-3073" "TMe-3267" "TMe-380"  "TMe-3780" "TMe-1579" "TMe-107" 
#>  [7] "TMe-108"  "TMe-1364" "TMe-207"  "TMe-1027" "TMe-387"  "TMe-2954"
#> [13] "TMe-82"   "TMe-3541" "TMe-386"  "TMe-3212" "TMe-396"  "TMe-1726"
#> [19] "TMe-3435" "TMe-3168" "TMe-513"  "TMe-840"  "TMe-1139" "TMe-3019"
#> [25] "TMe-802"  "TMe-1118" "TMe-2399" "TMe-3585" "TMe-885"  "TMe-1903"
#> [31] "TMe-2928" "TMe-3068" "TMe-3253" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-2643" "TMe-307"  "TMe-1127" "TMe-2835" "TMe-212"  "TMe-1011"
#>  [7] "TMe-813"  "TMe-305"  "TMe-1160" "TMe-1399" "TMe-1879" "TMe-1358"
#> [13] "TMe-826"  "TMe-474"  "TMe-982"  "TMe-747"  "TMe-1099" "TMe-629" 
#> [19] "TMe-417"  "TMe-685"  "TMe-217"  "TMe-667"  "TMe-2032" "TMe-868" 
#> [25] "TMe-2845" "TMe-1234" "TMe-730"  "TMe-341"  "TMe-2355" "TMe-2290"
#> [31] "TMe-1188" "TMe-2933" "TMe-1901" "TMe-2753" "TMe-1733" "TMe-850" 
#> [37] "TMe-2530" "TMe-2435" "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-1985" "TMe-1017" "TMe-1217" "TMe-1512" "TMe-1403"
#>  [7] "TMe-2308" "TMe-838"  "TMe-3708" "TMe-281"  "TMe-2535" "TMe-1164"
#> [13] "TMe-1216" "TMe-2963" "TMe-625"  "TMe-1661" "TMe-2389"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out8,
          use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "ward.D2")
```
