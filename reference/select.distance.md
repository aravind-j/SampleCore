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
counts <- c(I = 61, II = 41, III = 37, IV = 81, V = 80, VI = 37)

mand_accns <-
  c("TMe-34", "TMe-3423", "TMe-2018", "TMe-801", "TMe-551")

gp_vec <- setNames(as.character(data[, "Cluster"]), data[, "genotypes"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by centrality based methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Medoid-like Representative Sampling by Minimal Mean Distance
sel_mean_medoid_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "mean.medoid")
sel_mean_medoid_out
#> $I
#>  [1] "TMe-1696" "TMe-2251" "TMe-3003" "TMe-3202" "TMe-2388" "TMe-1226"
#>  [7] "TMe-1823" "TMe-3416" "TMe-2255" "TMe-990"  "TMe-677"  "TMe-2152"
#> [13] "TMe-1717" "TMe-1940" "TMe-2083" "TMe-1914" "TMe-898"  "TMe-2066"
#> [19] "TMe-3145" "TMe-1448" "TMe-3490" "TMe-1344" "TMe-486"  "TMe-2949"
#> [25] "TMe-3127" "TMe-716"  "TMe-933"  "TMe-2964" "TMe-3001" "TMe-3477"
#> [31] "TMe-3104" "TMe-3577" "TMe-2976" "TMe-20"   "TMe-1955" "TMe-3325"
#> [37] "TMe-1906" "TMe-1886" "TMe-3501" "TMe-882"  "TMe-3031" "TMe-3142"
#> [43] "TMe-1451" "TMe-2936" "TMe-3027" "TMe-910"  "TMe-1341" "TMe-489" 
#> [49] "TMe-469"  "TMe-3425" "TMe-3441" "TMe-3074" "TMe-501"  "TMe-1533"
#> [55] "TMe-1568" "TMe-2453" "TMe-896"  "TMe-1930" "TMe-2103" "TMe-2191"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-433"  "TMe-1969" "TMe-3203" "TMe-3365" "TMe-2242"
#>  [7] "TMe-3617" "TMe-1833" "TMe-377"  "TMe-3540" "TMe-455"  "TMe-2715"
#> [13] "TMe-1323" "TMe-2935" "TMe-1187" "TMe-3671" "TMe-3239" "TMe-2950"
#> [19] "TMe-3146" "TMe-1251" "TMe-3187" "TMe-2765" "TMe-1409" "TMe-2951"
#> [25] "TMe-1766" "TMe-1271" "TMe-1795" "TMe-2705" "TMe-1172" "TMe-2814"
#> [31] "TMe-3642" "TMe-1461" "TMe-2611" "TMe-1279" "TMe-890"  "TMe-3308"
#> [37] "TMe-2048" "TMe-3443" "TMe-1474" "TMe-2428" "TMe-1698"
#> 
#> $III
#>  [1] "TMe-2154" "TMe-2285" "TMe-3113" "TMe-2169" "TMe-3069" "TMe-3118"
#>  [7] "TMe-304"  "TMe-1138" "TMe-1725" "TMe-3326" "TMe-1836" "TMe-237" 
#> [13] "TMe-1797" "TMe-617"  "TMe-3324" "TMe-3028" "TMe-1675" "TMe-1792"
#> [19] "TMe-3317" "TMe-2748" "TMe-3216" "TMe-203"  "TMe-3731" "TMe-1421"
#> [25] "TMe-3085" "TMe-3166" "TMe-3005" "TMe-3207" "TMe-1965" "TMe-3335"
#> [31] "TMe-3147" "TMe-3592" "TMe-2356" "TMe-3010" "TMe-3796" "TMe-1262"
#> [37] "TMe-14"  
#> 
#> $IV
#>  [1] "TMe-1726" "TMe-2340" "TMe-3231" "TMe-1652" "TMe-3161" "TMe-2332"
#>  [7] "TMe-3327" "TMe-1867" "TMe-1211" "TMe-2363" "TMe-3763" "TMe-403" 
#> [13] "TMe-3233" "TMe-1828" "TMe-3562" "TMe-1167" "TMe-1651" "TMe-1553"
#> [19] "TMe-3218" "TMe-3168" "TMe-2337" "TMe-3412" "TMe-1118" "TMe-3212"
#> [25] "TMe-2210" "TMe-1948" "TMe-3214" "TMe-3494" "TMe-1996" "TMe-3660"
#> [31] "TMe-2390" "TMe-3232" "TMe-380"  "TMe-1765" "TMe-3108" "TMe-1166"
#> [37] "TMe-3107" "TMe-840"  "TMe-207"  "TMe-891"  "TMe-65"   "TMe-1402"
#> [43] "TMe-641"  "TMe-734"  "TMe-519"  "TMe-428"  "TMe-1369" "TMe-87"  
#> [49] "TMe-1959" "TMe-3073" "TMe-81"   "TMe-2946" "TMe-3780" "TMe-1470"
#> [55] "TMe-1325" "TMe-286"  "TMe-3378" "TMe-1128" "TMe-1139" "TMe-3779"
#> [61] "TMe-1374" "TMe-62"   "TMe-3602" "TMe-2367" "TMe-2399" "TMe-525" 
#> [67] "TMe-2988" "TMe-54"   "TMe-802"  "TMe-57"   "TMe-3535" "TMe-2989"
#> [73] "TMe-1246" "TMe-204"  "TMe-885"  "TMe-2378" "TMe-3290" "TMe-2039"
#> [79] "TMe-1397" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-456"  "TMe-474"  "TMe-1414" "TMe-574"  "TMe-2400"
#>  [7] "TMe-2845" "TMe-340"  "TMe-363"  "TMe-1268" "TMe-1160" "TMe-945" 
#> [13] "TMe-565"  "TMe-1234" "TMe-312"  "TMe-1979" "TMe-417"  "TMe-473" 
#> [19] "TMe-359"  "TMe-472"  "TMe-685"  "TMe-647"  "TMe-1862" "TMe-211" 
#> [25] "TMe-424"  "TMe-2271" "TMe-1183" "TMe-1622" "TMe-479"  "TMe-645" 
#> [31] "TMe-245"  "TMe-168"  "TMe-430"  "TMe-1559" "TMe-939"  "TMe-1760"
#> [37] "TMe-747"  "TMe-1382" "TMe-771"  "TMe-1715" "TMe-2361" "TMe-360" 
#> [43] "TMe-668"  "TMe-1299" "TMe-745"  "TMe-800"  "TMe-2009" "TMe-137" 
#> [49] "TMe-868"  "TMe-1572" "TMe-835"  "TMe-932"  "TMe-669"  "TMe-307" 
#> [55] "TMe-1273" "TMe-1345" "TMe-2441" "TMe-382"  "TMe-924"  "TMe-1308"
#> [61] "TMe-1100" "TMe-797"  "TMe-1521" "TMe-406"  "TMe-822"  "TMe-333" 
#> [67] "TMe-142"  "TMe-330"  "TMe-1188" "TMe-618"  "TMe-1609" "TMe-543" 
#> [73] "TMe-2057" "TMe-1522" "TMe-629"  "TMe-2439" "TMe-2124" "TMe-870" 
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-662"  "TMe-791"  "TMe-1847" "TMe-1164" "TMe-907" 
#>  [7] "TMe-1704" "TMe-222"  "TMe-1256" "TMe-985"  "TMe-726"  "TMe-1548"
#> [13] "TMe-1074" "TMe-1481" "TMe-809"  "TMe-720"  "TMe-752"  "TMe-1477"
#> [19] "TMe-1110" "TMe-273"  "TMe-2572" "TMe-936"  "TMe-1945" "TMe-683" 
#> [25] "TMe-1721" "TMe-652"  "TMe-2402" "TMe-505"  "TMe-1442" "TMe-1676"
#> [31] "TMe-2389" "TMe-781"  "TMe-1614" "TMe-678"  "TMe-3708" "TMe-814" 
#> [37] "TMe-281" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_mean_medoid_out, use.names = FALSE)) +
  labs(title = "mean.medoid")


# Medoid-like Representative Sampling by Minimal Median Distance
sel_median_medoid_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "median.medoid")
sel_median_medoid_out
#> $I
#>  [1] "TMe-3202" "TMe-1696" "TMe-2251" "TMe-3003" "TMe-677"  "TMe-1344"
#>  [7] "TMe-1823" "TMe-2255" "TMe-1226" "TMe-2388" "TMe-3416" "TMe-1451"
#> [13] "TMe-898"  "TMe-3127" "TMe-486"  "TMe-990"  "TMe-1940" "TMe-3145"
#> [19] "TMe-2066" "TMe-2949" "TMe-3001" "TMe-1914" "TMe-1717" "TMe-1955"
#> [25] "TMe-910"  "TMe-2083" "TMe-1448" "TMe-2152" "TMe-716"  "TMe-3142"
#> [31] "TMe-3477" "TMe-469"  "TMe-1906" "TMe-882"  "TMe-896"  "TMe-1341"
#> [37] "TMe-1533" "TMe-3104" "TMe-1886" "TMe-2976" "TMe-3499" "TMe-3490"
#> [43] "TMe-2939" "TMe-1568" "TMe-3534" "TMe-2191" "TMe-20"   "TMe-933" 
#> [49] "TMe-1777" "TMe-3031" "TMe-489"  "TMe-3460" "TMe-3623" "TMe-3074"
#> [55] "TMe-3170" "TMe-501"  "TMe-3501" "TMe-3548" "TMe-3027" "TMe-3577"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3365" "TMe-3495" "TMe-3203" "TMe-3540" "TMe-1833" "TMe-2935"
#>  [7] "TMe-1795" "TMe-2715" "TMe-377"  "TMe-455"  "TMe-433"  "TMe-1766"
#> [13] "TMe-2705" "TMe-2242" "TMe-1271" "TMe-1969" "TMe-2814" "TMe-2765"
#> [19] "TMe-3671" "TMe-1323" "TMe-1187" "TMe-3187" "TMe-2950" "TMe-1172"
#> [25] "TMe-1409" "TMe-2611" "TMe-3239" "TMe-2951" "TMe-2414" "TMe-1505"
#> [31] "TMe-1251" "TMe-1279" "TMe-3617" "TMe-3146" "TMe-3498" "TMe-1461"
#> [37] "TMe-3443" "TMe-251"  "TMe-3020" "TMe-2048" "TMe-3612"
#> 
#> $III
#>  [1] "TMe-2285" "TMe-3118" "TMe-2154" "TMe-3069" "TMe-1836" "TMe-3113"
#>  [7] "TMe-2169" "TMe-3028" "TMe-1138" "TMe-304"  "TMe-1725" "TMe-3317"
#> [13] "TMe-3324" "TMe-237"  "TMe-1797" "TMe-1421" "TMe-3731" "TMe-3085"
#> [19] "TMe-3796" "TMe-3326" "TMe-3216" "TMe-3147" "TMe-617"  "TMe-203" 
#> [25] "TMe-1792" "TMe-3207" "TMe-2356" "TMe-2748" "TMe-3166" "TMe-1675"
#> [31] "TMe-3592" "TMe-1956" "TMe-64"   "TMe-1910" "TMe-14"   "TMe-3005"
#> [37] "TMe-3321"
#> 
#> $IV
#>  [1] "TMe-3161" "TMe-3231" "TMe-1726" "TMe-1652" "TMe-2340" "TMe-3562"
#>  [7] "TMe-3212" "TMe-1211" "TMe-1828" "TMe-3327" "TMe-2363" "TMe-1167"
#> [13] "TMe-3660" "TMe-1867" "TMe-3494" "TMe-2337" "TMe-1996" "TMe-1948"
#> [19] "TMe-3412" "TMe-3763" "TMe-2332" "TMe-403"  "TMe-3233" "TMe-3602"
#> [25] "TMe-1553" "TMe-380"  "TMe-734"  "TMe-3218" "TMe-2210" "TMe-1651"
#> [31] "TMe-1118" "TMe-1765" "TMe-1325" "TMe-2390" "TMe-3779" "TMe-1166"
#> [37] "TMe-3073" "TMe-87"   "TMe-3108" "TMe-65"   "TMe-81"   "TMe-2946"
#> [43] "TMe-428"  "TMe-3232" "TMe-3780" "TMe-3290" "TMe-3168" "TMe-3107"
#> [49] "TMe-2367" "TMe-1369" "TMe-891"  "TMe-641"  "TMe-525"  "TMe-840" 
#> [55] "TMe-207"  "TMe-1402" "TMe-396"  "TMe-3214" "TMe-1959" "TMe-3535"
#> [61] "TMe-3206" "TMe-1128" "TMe-3191" "TMe-1903" "TMe-3757" "TMe-204" 
#> [67] "TMe-286"  "TMe-519"  "TMe-30"   "TMe-2375" "TMe-2410" "TMe-3538"
#> [73] "TMe-3257" "TMe-2988" "TMe-1139" "TMe-54"   "TMe-57"   "TMe-1246"
#> [79] "TMe-1374" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-456"  "TMe-474"  "TMe-574"  "TMe-2845" "TMe-1414"
#>  [7] "TMe-2400" "TMe-565"  "TMe-340"  "TMe-1268" "TMe-312"  "TMe-417" 
#> [13] "TMe-472"  "TMe-363"  "TMe-1160" "TMe-359"  "TMe-945"  "TMe-1979"
#> [19] "TMe-1234" "TMe-647"  "TMe-685"  "TMe-211"  "TMe-747"  "TMe-473" 
#> [25] "TMe-2271" "TMe-745"  "TMe-645"  "TMe-1183" "TMe-1622" "TMe-1760"
#> [31] "TMe-1559" "TMe-1862" "TMe-1308" "TMe-245"  "TMe-479"  "TMe-2361"
#> [37] "TMe-1273" "TMe-1715" "TMe-939"  "TMe-430"  "TMe-424"  "TMe-1382"
#> [43] "TMe-800"  "TMe-771"  "TMe-360"  "TMe-168"  "TMe-2009" "TMe-1572"
#> [49] "TMe-868"  "TMe-668"  "TMe-137"  "TMe-1345" "TMe-1299" "TMe-618" 
#> [55] "TMe-2441" "TMe-1609" "TMe-307"  "TMe-932"  "TMe-382"  "TMe-669" 
#> [61] "TMe-797"  "TMe-333"  "TMe-1880" "TMe-1188" "TMe-1957" "TMe-822" 
#> [67] "TMe-323"  "TMe-835"  "TMe-1100" "TMe-543"  "TMe-142"  "TMe-432" 
#> [73] "TMe-977"  "TMe-330"  "TMe-927"  "TMe-1227" "TMe-629"  "TMe-2124"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-907"  "TMe-1847" "TMe-1164" "TMe-662"  "TMe-1704"
#>  [7] "TMe-791"  "TMe-1256" "TMe-809"  "TMe-1110" "TMe-222"  "TMe-985" 
#> [13] "TMe-1548" "TMe-1945" "TMe-726"  "TMe-1074" "TMe-683"  "TMe-1477"
#> [19] "TMe-781"  "TMe-2572" "TMe-936"  "TMe-1721" "TMe-1481" "TMe-652" 
#> [25] "TMe-1676" "TMe-720"  "TMe-752"  "TMe-1217" "TMe-2389" "TMe-1442"
#> [31] "TMe-814"  "TMe-505"  "TMe-273"  "TMe-2402" "TMe-856"  "TMe-678" 
#> [37] "TMe-1614"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_median_medoid_out, use.names = FALSE)) +
  labs(title = "median.medoid")


# Representative Sampling by Proximity to Group Centroid
sel_group_centroid_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "nearest.centroid")
sel_group_centroid_out
#> $I
#>  [1] "TMe-1696" "TMe-2251" "TMe-2964" "TMe-3202" "TMe-2388" "TMe-1226"
#>  [7] "TMe-1823" "TMe-3416" "TMe-2255" "TMe-677"  "TMe-990"  "TMe-2152"
#> [13] "TMe-1717" "TMe-2083" "TMe-1940" "TMe-1914" "TMe-898"  "TMe-2066"
#> [19] "TMe-3145" "TMe-1448" "TMe-3490" "TMe-486"  "TMe-1344" "TMe-2949"
#> [25] "TMe-3127" "TMe-2943" "TMe-716"  "TMe-933"  "TMe-3001" "TMe-3477"
#> [31] "TMe-3104" "TMe-2976" "TMe-3548" "TMe-3729" "TMe-1955" "TMe-1906"
#> [37] "TMe-3325" "TMe-3493" "TMe-882"  "TMe-1886" "TMe-3142" "TMe-3031"
#> [43] "TMe-1451" "TMe-2936" "TMe-3027" "TMe-33"   "TMe-1341" "TMe-489" 
#> [49] "TMe-469"  "TMe-3424" "TMe-3437" "TMe-501"  "TMe-3074" "TMe-1533"
#> [55] "TMe-1568" "TMe-2453" "TMe-896"  "TMe-1930" "TMe-2103" "TMe-2191"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-433"  "TMe-1969" "TMe-3203" "TMe-3365" "TMe-2242"
#>  [7] "TMe-3617" "TMe-1833" "TMe-377"  "TMe-3540" "TMe-455"  "TMe-2715"
#> [13] "TMe-1323" "TMe-2935" "TMe-1187" "TMe-3671" "TMe-3239" "TMe-2950"
#> [19] "TMe-3146" "TMe-1251" "TMe-3187" "TMe-2765" "TMe-1409" "TMe-2951"
#> [25] "TMe-1766" "TMe-1271" "TMe-1795" "TMe-2705" "TMe-1172" "TMe-2814"
#> [31] "TMe-3642" "TMe-1461" "TMe-2611" "TMe-1279" "TMe-890"  "TMe-3308"
#> [37] "TMe-2048" "TMe-3443" "TMe-1474" "TMe-2428" "TMe-1698"
#> 
#> $III
#>  [1] "TMe-2154" "TMe-2285" "TMe-3113" "TMe-2169" "TMe-3069" "TMe-3118"
#>  [7] "TMe-304"  "TMe-1138" "TMe-1725" "TMe-3326" "TMe-1836" "TMe-237" 
#> [13] "TMe-1797" "TMe-617"  "TMe-3324" "TMe-3028" "TMe-1675" "TMe-1792"
#> [19] "TMe-3317" "TMe-2748" "TMe-3216" "TMe-203"  "TMe-3731" "TMe-1421"
#> [25] "TMe-3085" "TMe-3166" "TMe-3005" "TMe-3207" "TMe-1965" "TMe-3335"
#> [31] "TMe-3147" "TMe-3592" "TMe-2356" "TMe-3010" "TMe-3796" "TMe-1262"
#> [37] "TMe-14"  
#> 
#> $IV
#>  [1] "TMe-1700" "TMe-2332" "TMe-1597" "TMe-3144" "TMe-3218" "TMe-2036"
#>  [7] "TMe-3297" "TMe-1820" "TMe-1179" "TMe-2247" "TMe-3760" "TMe-403" 
#> [13] "TMe-3231" "TMe-3541" "TMe-1652" "TMe-1166" "TMe-3206" "TMe-1580"
#> [19] "TMe-1526" "TMe-3161" "TMe-2318" "TMe-1106" "TMe-2172" "TMe-3390"
#> [25] "TMe-3198" "TMe-1923" "TMe-3451" "TMe-3044" "TMe-3225" "TMe-3619"
#> [31] "TMe-1987" "TMe-1651" "TMe-380"  "TMe-2375" "TMe-3097" "TMe-1162"
#> [37] "TMe-3045" "TMe-834"  "TMe-1376" "TMe-641"  "TMe-190"  "TMe-519" 
#> [43] "TMe-840"  "TMe-734"  "TMe-54"   "TMe-428"  "TMe-1368" "TMe-81"  
#> [49] "TMe-3072" "TMe-3761" "TMe-1948" "TMe-1419" "TMe-2924" "TMe-72"  
#> [55] "TMe-1229" "TMe-280"  "TMe-3327" "TMe-1123" "TMe-1094" "TMe-1369"
#> [61] "TMe-3773" "TMe-62"   "TMe-2340" "TMe-2378" "TMe-3538" "TMe-2981"
#> [67] "TMe-525"  "TMe-761"  "TMe-38"   "TMe-57"   "TMe-204"  "TMe-3511"
#> [73] "TMe-802"  "TMe-1231" "TMe-2946" "TMe-1336" "TMe-2025" "TMe-3276"
#> [79] "TMe-2367" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-456"  "TMe-474"  "TMe-1411" "TMe-565"  "TMe-2355"
#>  [7] "TMe-2769" "TMe-340"  "TMe-363"  "TMe-1257" "TMe-1127" "TMe-939" 
#> [13] "TMe-551"  "TMe-1220" "TMe-287"  "TMe-1963" "TMe-417"  "TMe-473" 
#> [19] "TMe-472"  "TMe-359"  "TMe-669"  "TMe-197"  "TMe-1859" "TMe-623" 
#> [25] "TMe-2003" "TMe-424"  "TMe-1534" "TMe-1160" "TMe-479"  "TMe-217" 
#> [31] "TMe-629"  "TMe-142"  "TMe-1551" "TMe-370"  "TMe-932"  "TMe-1730"
#> [37] "TMe-731"  "TMe-2339" "TMe-1381" "TMe-764"  "TMe-1694" "TMe-360" 
#> [43] "TMe-667"  "TMe-1293" "TMe-1979" "TMe-786"  "TMe-682"  "TMe-2973"
#> [49] "TMe-861"  "TMe-1556" "TMe-826"  "TMe-929"  "TMe-307"  "TMe-655" 
#> [55] "TMe-1265" "TMe-1295" "TMe-2435" "TMe-382"  "TMe-889"  "TMe-1305"
#> [61] "TMe-1099" "TMe-1500" "TMe-795"  "TMe-406"  "TMe-326"  "TMe-3185"
#> [67] "TMe-813"  "TMe-323"  "TMe-1158" "TMe-1559" "TMe-2050" "TMe-543" 
#> [73] "TMe-612"  "TMe-1507" "TMe-603"  "TMe-2425" "TMe-1629" "TMe-2119"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-662"  "TMe-791"  "TMe-1847" "TMe-1164" "TMe-907" 
#>  [7] "TMe-1704" "TMe-222"  "TMe-1256" "TMe-985"  "TMe-726"  "TMe-1548"
#> [13] "TMe-1074" "TMe-1481" "TMe-809"  "TMe-720"  "TMe-752"  "TMe-1477"
#> [19] "TMe-1110" "TMe-273"  "TMe-2572" "TMe-936"  "TMe-1945" "TMe-683" 
#> [25] "TMe-1721" "TMe-652"  "TMe-2402" "TMe-505"  "TMe-1442" "TMe-1676"
#> [31] "TMe-2389" "TMe-781"  "TMe-1614" "TMe-678"  "TMe-3708" "TMe-814" 
#> [37] "TMe-281" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_group_centroid_out, use.names = FALSE)) +
  labs(title = "nearest.centroid")


# Representative Sampling by Proximity to Group Spatial Median
sel_group_median_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "nearest.median")
sel_group_median_out
#> $I
#>  [1] "TMe-1696" "TMe-2251" "TMe-2964" "TMe-3202" "TMe-2388" "TMe-1226"
#>  [7] "TMe-1823" "TMe-3416" "TMe-2255" "TMe-2152" "TMe-677"  "TMe-1717"
#> [13] "TMe-990"  "TMe-2083" "TMe-1940" "TMe-2066" "TMe-898"  "TMe-3145"
#> [19] "TMe-1914" "TMe-1448" "TMe-486"  "TMe-1344" "TMe-3127" "TMe-933" 
#> [25] "TMe-3490" "TMe-716"  "TMe-2949" "TMe-3001" "TMe-2943" "TMe-1955"
#> [31] "TMe-3729" "TMe-2976" "TMe-1451" "TMe-3477" "TMe-1906" "TMe-3548"
#> [37] "TMe-3325" "TMe-3142" "TMe-3031" "TMe-3104" "TMe-3027" "TMe-882" 
#> [43] "TMe-33"   "TMe-3493" "TMe-1886" "TMe-501"  "TMe-469"  "TMe-1341"
#> [49] "TMe-489"  "TMe-2936" "TMe-1533" "TMe-896"  "TMe-2453" "TMe-1568"
#> [55] "TMe-3437" "TMe-3074" "TMe-1930" "TMe-3424" "TMe-2103" "TMe-2191"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-433"  "TMe-1969" "TMe-3203" "TMe-2242" "TMe-3365"
#>  [7] "TMe-1833" "TMe-3540" "TMe-3617" "TMe-377"  "TMe-455"  "TMe-2715"
#> [13] "TMe-1323" "TMe-2935" "TMe-3671" "TMe-1187" "TMe-3239" "TMe-3146"
#> [19] "TMe-1251" "TMe-2950" "TMe-2765" "TMe-1409" "TMe-3187" "TMe-1766"
#> [25] "TMe-1795" "TMe-2951" "TMe-2705" "TMe-1271" "TMe-2814" "TMe-1172"
#> [31] "TMe-1461" "TMe-2611" "TMe-3642" "TMe-1279" "TMe-890"  "TMe-3308"
#> [37] "TMe-2048" "TMe-1474" "TMe-3443" "TMe-2428" "TMe-681" 
#> 
#> $III
#>  [1] "TMe-2154" "TMe-3113" "TMe-2285" "TMe-2169" "TMe-3069" "TMe-3118"
#>  [7] "TMe-1725" "TMe-304"  "TMe-1138" "TMe-3326" "TMe-1836" "TMe-237" 
#> [13] "TMe-1797" "TMe-3324" "TMe-617"  "TMe-3028" "TMe-1792" "TMe-1675"
#> [19] "TMe-3317" "TMe-3216" "TMe-2748" "TMe-1421" "TMe-203"  "TMe-3731"
#> [25] "TMe-3085" "TMe-3166" "TMe-3005" "TMe-3335" "TMe-1965" "TMe-3207"
#> [31] "TMe-3147" "TMe-3010" "TMe-2356" "TMe-3592" "TMe-3796" "TMe-14"  
#> [37] "TMe-1262"
#> 
#> $IV
#>  [1] "TMe-1700" "TMe-2332" "TMe-1597" "TMe-3144" "TMe-3218" "TMe-2036"
#>  [7] "TMe-3297" "TMe-1820" "TMe-1179" "TMe-2247" "TMe-403"  "TMe-3760"
#> [13] "TMe-3231" "TMe-3541" "TMe-1652" "TMe-1580" "TMe-3206" "TMe-1526"
#> [19] "TMe-1166" "TMe-2318" "TMe-1106" "TMe-2172" "TMe-3390" "TMe-1923"
#> [25] "TMe-3198" "TMe-3161" "TMe-1651" "TMe-3451" "TMe-3225" "TMe-380" 
#> [31] "TMe-3619" "TMe-3044" "TMe-1987" "TMe-2375" "TMe-3097" "TMe-1162"
#> [37] "TMe-834"  "TMe-3045" "TMe-1376" "TMe-734"  "TMe-641"  "TMe-54"  
#> [43] "TMe-519"  "TMe-190"  "TMe-840"  "TMe-428"  "TMe-81"   "TMe-3072"
#> [49] "TMe-2924" "TMe-1419" "TMe-1368" "TMe-3761" "TMe-1948" "TMe-72"  
#> [55] "TMe-1229" "TMe-280"  "TMe-3327" "TMe-1369" "TMe-1094" "TMe-3773"
#> [61] "TMe-1123" "TMe-2340" "TMe-62"   "TMe-38"   "TMe-57"   "TMe-2981"
#> [67] "TMe-3538" "TMe-3511" "TMe-204"  "TMe-2378" "TMe-761"  "TMe-802" 
#> [73] "TMe-1336" "TMe-525"  "TMe-2946" "TMe-1231" "TMe-27"   "TMe-3276"
#> [79] "TMe-3055" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-332"  "TMe-456"  "TMe-474"  "TMe-1411" "TMe-565"  "TMe-2769"
#>  [7] "TMe-2355" "TMe-340"  "TMe-363"  "TMe-1127" "TMe-1257" "TMe-939" 
#> [13] "TMe-551"  "TMe-287"  "TMe-1220" "TMe-1963" "TMe-473"  "TMe-417" 
#> [19] "TMe-359"  "TMe-472"  "TMe-197"  "TMe-669"  "TMe-1859" "TMe-1534"
#> [25] "TMe-623"  "TMe-1160" "TMe-2003" "TMe-424"  "TMe-629"  "TMe-217" 
#> [31] "TMe-479"  "TMe-2339" "TMe-370"  "TMe-932"  "TMe-1551" "TMe-731" 
#> [37] "TMe-764"  "TMe-1694" "TMe-360"  "TMe-142"  "TMe-1730" "TMe-1381"
#> [43] "TMe-682"  "TMe-1293" "TMe-667"  "TMe-786"  "TMe-1556" "TMe-1979"
#> [49] "TMe-2973" "TMe-861"  "TMe-826"  "TMe-1265" "TMe-2435" "TMe-655" 
#> [55] "TMe-307"  "TMe-1295" "TMe-929"  "TMe-382"  "TMe-1500" "TMe-795" 
#> [61] "TMe-889"  "TMe-1099" "TMe-1305" "TMe-326"  "TMe-323"  "TMe-3185"
#> [67] "TMe-813"  "TMe-2425" "TMe-406"  "TMe-543"  "TMe-1559" "TMe-612" 
#> [73] "TMe-2050" "TMe-863"  "TMe-603"  "TMe-1158" "TMe-1322" "TMe-1507"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1550" "TMe-662"  "TMe-791"  "TMe-1847" "TMe-1164" "TMe-907" 
#>  [7] "TMe-1704" "TMe-222"  "TMe-1256" "TMe-985"  "TMe-726"  "TMe-1074"
#> [13] "TMe-1548" "TMe-1481" "TMe-809"  "TMe-720"  "TMe-1110" "TMe-752" 
#> [19] "TMe-1477" "TMe-273"  "TMe-1945" "TMe-936"  "TMe-1721" "TMe-2572"
#> [25] "TMe-683"  "TMe-505"  "TMe-2402" "TMe-1442" "TMe-781"  "TMe-652" 
#> [31] "TMe-1676" "TMe-2389" "TMe-1614" "TMe-678"  "TMe-3708" "TMe-814" 
#> [37] "TMe-281" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_group_median_out, use.names = FALSE)) +
  labs(title = "nearest.median")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by peripheral/extremity based methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Peripheral Sampling by Maximal Mean Distance
sel_mean_peripheral_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "mean.peripheral")
sel_mean_peripheral_out
#> $I
#>  [1] "TMe-2967" "TMe-3685" "TMe-3705" "TMe-3292" "TMe-3223" "TMe-3667"
#>  [7] "TMe-41"   "TMe-3475" "TMe-606"  "TMe-1425" "TMe-3282" "TMe-3319"
#> [13] "TMe-3392" "TMe-2965" "TMe-3437" "TMe-1564" "TMe-2943" "TMe-2955"
#> [19] "TMe-2993" "TMe-756"  "TMe-2069" "TMe-717"  "TMe-3415" "TMe-3396"
#> [25] "TMe-3389" "TMe-500"  "TMe-3065" "TMe-2996" "TMe-867"  "TMe-3641"
#> [31] "TMe-2779" "TMe-3115" "TMe-3163" "TMe-841"  "TMe-3151" "TMe-3433"
#> [37] "TMe-3249" "TMe-1360" "TMe-1190" "TMe-3330" "TMe-410"  "TMe-3266"
#> [43] "TMe-2027" "TMe-3102" "TMe-3398" "TMe-3264" "TMe-2604" "TMe-815" 
#> [49] "TMe-2785" "TMe-3478" "TMe-1218" "TMe-727"  "TMe-1333" "TMe-3323"
#> [55] "TMe-569"  "TMe-3471" "TMe-2985" "TMe-2984" "TMe-3281" "TMe-1466"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3466" "TMe-3338" "TMe-47"   "TMe-2860" "TMe-1919" "TMe-3605"
#>  [7] "TMe-1907" "TMe-86"   "TMe-1107" "TMe-2128" "TMe-369"  "TMe-2033"
#> [13] "TMe-2952" "TMe-2056" "TMe-3141" "TMe-3805" "TMe-2568" "TMe-2329"
#> [19] "TMe-1385" "TMe-1137" "TMe-3766" "TMe-3800" "TMe-674"  "TMe-3210"
#> [25] "TMe-3200" "TMe-2352" "TMe-2891" "TMe-3690" "TMe-3140" "TMe-509" 
#> [31] "TMe-1619" "TMe-1444" "TMe-478"  "TMe-404"  "TMe-2866" "TMe-171" 
#> [37] "TMe-289"  "TMe-2204" "TMe-160"  "TMe-477"  "TMe-3222"
#> 
#> $III
#>  [1] "TMe-3596" "TMe-2926" "TMe-1897" "TMe-3274" "TMe-3565" "TMe-267" 
#>  [7] "TMe-3804" "TMe-3299" "TMe-785"  "TMe-234"  "TMe-635"  "TMe-261" 
#> [13] "TMe-1868" "TMe-3445" "TMe-35"   "TMe-3663" "TMe-116"  "TMe-773" 
#> [19] "TMe-425"  "TMe-3631" "TMe-13"   "TMe-381"  "TMe-70"   "TMe-187" 
#> [25] "TMe-3032" "TMe-2823" "TMe-3383" "TMe-2897" "TMe-4"    "TMe-141" 
#> [31] "TMe-3143" "TMe-3043" "TMe-3715" "TMe-946"  "TMe-3048" "TMe-2756"
#> [37] "TMe-2811"
#> 
#> $IV
#>  [1] "TMe-2809" "TMe-241"  "TMe-1434" "TMe-550"  "TMe-812"  "TMe-761" 
#>  [7] "TMe-3417" "TMe-2924" "TMe-3442" "TMe-1988" "TMe-1525" "TMe-698" 
#> [13] "TMe-1020" "TMe-3730" "TMe-1078" "TMe-3054" "TMe-608"  "TMe-3273"
#> [19] "TMe-3573" "TMe-2043" "TMe-226"  "TMe-180"  "TMe-2956" "TMe-3275"
#> [25] "TMe-2567" "TMe-154"  "TMe-3545" "TMe-421"  "TMe-1229" "TMe-2971"
#> [31] "TMe-317"  "TMe-3025" "TMe-3707" "TMe-650"  "TMe-3167" "TMe-1162"
#> [37] "TMe-3189" "TMe-3255" "TMe-1129" "TMe-3297" "TMe-27"   "TMe-427" 
#> [43] "TMe-1419" "TMe-3256" "TMe-3055" "TMe-3089" "TMe-2958" "TMe-766" 
#> [49] "TMe-388"  "TMe-372"  "TMe-737"  "TMe-2979" "TMe-225"  "TMe-3701"
#> [55] "TMe-684"  "TMe-280"  "TMe-2552" "TMe-1376" "TMe-2240" "TMe-3225"
#> [61] "TMe-1328" "TMe-209"  "TMe-1479" "TMe-3053" "TMe-1377" "TMe-259" 
#> [67] "TMe-3594" "TMe-1348" "TMe-3527" "TMe-2833" "TMe-3040" "TMe-619" 
#> [73] "TMe-824"  "TMe-656"  "TMe-1891" "TMe-962"  "TMe-540"  "TMe-184" 
#> [79] "TMe-787"  "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-3034" "TMe-2121" "TMe-3628" "TMe-419"  "TMe-3736" "TMe-3185"
#>  [7] "TMe-536"  "TMe-2853" "TMe-378"  "TMe-712"  "TMe-2355" "TMe-2589"
#> [13] "TMe-2050" "TMe-1458" "TMe-2531" "TMe-616"  "TMe-723"  "TMe-3356"
#> [19] "TMe-1381" "TMe-1099" "TMe-1455" "TMe-3332" "TMe-399"  "TMe-700" 
#> [25] "TMe-782"  "TMe-2612" "TMe-603"  "TMe-277"  "TMe-861"  "TMe-2953"
#> [31] "TMe-284"  "TMe-1391" "TMe-2907" "TMe-2820" "TMe-2296" "TMe-3411"
#> [37] "TMe-2413" "TMe-3408" "TMe-2933" "TMe-2790" "TMe-832"  "TMe-2973"
#> [43] "TMe-2435" "TMe-3363" "TMe-1257" "TMe-1399" "TMe-1730" "TMe-3311"
#> [49] "TMe-2003" "TMe-344"  "TMe-2905" "TMe-2060" "TMe-189"  "TMe-1042"
#> [55] "TMe-2290" "TMe-294"  "TMe-2688" "TMe-1924" "TMe-48"   "TMe-969" 
#> [61] "TMe-1357" "TMe-1159" "TMe-2750" "TMe-730"  "TMe-208"  "TMe-532" 
#> [67] "TMe-3659" "TMe-2139" "TMe-1311" "TMe-2761" "TMe-755"  "TMe-1963"
#> [73] "TMe-158"  "TMe-481"  "TMe-2863" "TMe-287"  "TMe-2904" "TMe-2307"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-3095" "TMe-3549" "TMe-2983" "TMe-2035" "TMe-3116" "TMe-1416"
#>  [7] "TMe-693"  "TMe-2963" "TMe-1232" "TMe-1769" "TMe-1403" "TMe-2064"
#> [13] "TMe-1875" "TMe-1124" "TMe-1383" "TMe-1518" "TMe-1428" "TMe-1174"
#> [19] "TMe-2510" "TMe-1775" "TMe-725"  "TMe-1509" "TMe-1992" "TMe-465" 
#> [25] "TMe-1035" "TMe-1985" "TMe-2791" "TMe-836"  "TMe-690"  "TMe-1261"
#> [31] "TMe-1182" "TMe-3387" "TMe-1506" "TMe-751"  "TMe-2534" "TMe-361" 
#> [37] "TMe-1053"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_mean_peripheral_out, use.names = FALSE)) +
  labs(title = "mean.peripheral")


# Peripheral Sampling by Maximal Median Distance
sel_median_peripheral_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "median.peripheral")
sel_median_peripheral_out
#> $I
#>  [1] "TMe-3685" "TMe-2967" "TMe-3705" "TMe-3223" "TMe-3292" "TMe-3667"
#>  [7] "TMe-3475" "TMe-41"   "TMe-3282" "TMe-3319" "TMe-1425" "TMe-606" 
#> [13] "TMe-3392" "TMe-2993" "TMe-2955" "TMe-2943" "TMe-2965" "TMe-1564"
#> [19] "TMe-756"  "TMe-3437" "TMe-2069" "TMe-3641" "TMe-3389" "TMe-3415"
#> [25] "TMe-3065" "TMe-3396" "TMe-717"  "TMe-2779" "TMe-3249" "TMe-3266"
#> [31] "TMe-867"  "TMe-3330" "TMe-2996" "TMe-3115" "TMe-841"  "TMe-500" 
#> [37] "TMe-2604" "TMe-3163" "TMe-3433" "TMe-1190" "TMe-3102" "TMe-3264"
#> [43] "TMe-1360" "TMe-3398" "TMe-3151" "TMe-410"  "TMe-2027" "TMe-3323"
#> [49] "TMe-1218" "TMe-569"  "TMe-727"  "TMe-2984" "TMe-3281" "TMe-715" 
#> [55] "TMe-815"  "TMe-3471" "TMe-3601" "TMe-2985" "TMe-3478" "TMe-2934"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3466" "TMe-47"   "TMe-3338" "TMe-1919" "TMe-2860" "TMe-3605"
#>  [7] "TMe-1907" "TMe-1107" "TMe-2128" "TMe-86"   "TMe-369"  "TMe-2056"
#> [13] "TMe-2952" "TMe-2033" "TMe-674"  "TMe-3141" "TMe-1385" "TMe-1137"
#> [19] "TMe-3805" "TMe-2568" "TMe-3690" "TMe-3800" "TMe-1619" "TMe-2329"
#> [25] "TMe-2891" "TMe-3766" "TMe-3210" "TMe-2352" "TMe-509"  "TMe-3200"
#> [31] "TMe-1242" "TMe-3140" "TMe-2866" "TMe-1444" "TMe-171"  "TMe-289" 
#> [37] "TMe-404"  "TMe-539"  "TMe-478"  "TMe-477"  "TMe-2204"
#> 
#> $III
#>  [1] "TMe-3596" "TMe-2926" "TMe-1897" "TMe-3274" "TMe-3565" "TMe-3299"
#>  [7] "TMe-267"  "TMe-116"  "TMe-234"  "TMe-261"  "TMe-785"  "TMe-3663"
#> [13] "TMe-3804" "TMe-773"  "TMe-3631" "TMe-35"   "TMe-635"  "TMe-1868"
#> [19] "TMe-3445" "TMe-425"  "TMe-13"   "TMe-3043" "TMe-70"   "TMe-2756"
#> [25] "TMe-3143" "TMe-4"    "TMe-2823" "TMe-187"  "TMe-381"  "TMe-2897"
#> [31] "TMe-1889" "TMe-3032" "TMe-3383" "TMe-3048" "TMe-3715" "TMe-141" 
#> [37] "TMe-1993"
#> 
#> $IV
#>  [1] "TMe-2809" "TMe-241"  "TMe-1434" "TMe-550"  "TMe-812"  "TMe-761" 
#>  [7] "TMe-2924" "TMe-3417" "TMe-1988" "TMe-1078" "TMe-3442" "TMe-1020"
#> [13] "TMe-608"  "TMe-698"  "TMe-1525" "TMe-3273" "TMe-3730" "TMe-3573"
#> [19] "TMe-3054" "TMe-317"  "TMe-226"  "TMe-2043" "TMe-2971" "TMe-2956"
#> [25] "TMe-2567" "TMe-1419" "TMe-1229" "TMe-154"  "TMe-180"  "TMe-3545"
#> [31] "TMe-650"  "TMe-421"  "TMe-3707" "TMe-388"  "TMe-1129" "TMe-1162"
#> [37] "TMe-3297" "TMe-3275" "TMe-427"  "TMe-3255" "TMe-3025" "TMe-3189"
#> [43] "TMe-3089" "TMe-27"   "TMe-372"  "TMe-684"  "TMe-3167" "TMe-3701"
#> [49] "TMe-3256" "TMe-2552" "TMe-2979" "TMe-766"  "TMe-1328" "TMe-2958"
#> [55] "TMe-3242" "TMe-3527" "TMe-1479" "TMe-737"  "TMe-1377" "TMe-1010"
#> [61] "TMe-1297" "TMe-787"  "TMe-3055" "TMe-2833" "TMe-280"  "TMe-1376"
#> [67] "TMe-619"  "TMe-3040" "TMe-225"  "TMe-259"  "TMe-962"  "TMe-656" 
#> [73] "TMe-353"  "TMe-1891" "TMe-3225" "TMe-824"  "TMe-2240" "TMe-209" 
#> [79] "TMe-3053" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-2121" "TMe-3034" "TMe-3628" "TMe-419"  "TMe-3185" "TMe-3736"
#>  [7] "TMe-536"  "TMe-1458" "TMe-378"  "TMe-2853" "TMe-1455" "TMe-1381"
#> [13] "TMe-2355" "TMe-2589" "TMe-712"  "TMe-2050" "TMe-723"  "TMe-616" 
#> [19] "TMe-3356" "TMe-3332" "TMe-2531" "TMe-1099" "TMe-700"  "TMe-1391"
#> [25] "TMe-2612" "TMe-782"  "TMe-277"  "TMe-399"  "TMe-2296" "TMe-2953"
#> [31] "TMe-2907" "TMe-603"  "TMe-2820" "TMe-861"  "TMe-284"  "TMe-2790"
#> [37] "TMe-832"  "TMe-2413" "TMe-3411" "TMe-2933" "TMe-3408" "TMe-3363"
#> [43] "TMe-1257" "TMe-2973" "TMe-189"  "TMe-1924" "TMe-2905" "TMe-2435"
#> [49] "TMe-1159" "TMe-344"  "TMe-1730" "TMe-2290" "TMe-1399" "TMe-48"  
#> [55] "TMe-3311" "TMe-1042" "TMe-294"  "TMe-2003" "TMe-2139" "TMe-969" 
#> [61] "TMe-2688" "TMe-730"  "TMe-532"  "TMe-287"  "TMe-2060" "TMe-1357"
#> [67] "TMe-755"  "TMe-2307" "TMe-158"  "TMe-1375" "TMe-3659" "TMe-1311"
#> [73] "TMe-208"  "TMe-481"  "TMe-2750" "TMe-1332" "TMe-2761" "TMe-2904"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-3095" "TMe-3549" "TMe-1403" "TMe-2983" "TMe-3116" "TMe-693" 
#>  [7] "TMe-1416" "TMe-2035" "TMe-1232" "TMe-1769" "TMe-2963" "TMe-2064"
#> [13] "TMe-1875" "TMe-1124" "TMe-1383" "TMe-1518" "TMe-1775" "TMe-2510"
#> [19] "TMe-1428" "TMe-1992" "TMe-1174" "TMe-1035" "TMe-1985" "TMe-725" 
#> [25] "TMe-465"  "TMe-1261" "TMe-2791" "TMe-1509" "TMe-690"  "TMe-751" 
#> [31] "TMe-836"  "TMe-1182" "TMe-361"  "TMe-2534" "TMe-2196" "TMe-269" 
#> [37] "TMe-130" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_median_peripheral_out, use.names = FALSE)) +
  labs(title = "median.peripheral")


# Peripheral Sampling by Maximal Eccentricity
sel_eccentricity_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "eccentricity")
sel_eccentricity_out
#> $I
#>  [1] "TMe-606"  "TMe-3667" "TMe-3030" "TMe-3705" "TMe-3685" "TMe-1533"
#>  [7] "TMe-3480" "TMe-3437" "TMe-2943" "TMe-132"  "TMe-3719" "TMe-1564"
#> [13] "TMe-3292" "TMe-41"   "TMe-2914" "TMe-3698" "TMe-2967" "TMe-3169"
#> [19] "TMe-3208" "TMe-2383" "TMe-1347" "TMe-2004" "TMe-815"  "TMe-3223"
#> [25] "TMe-2216" "TMe-896"  "TMe-1621" "TMe-3110" "TMe-2965" "TMe-1468"
#> [31] "TMe-1532" "TMe-1425" "TMe-2779" "TMe-3353" "TMe-1830" "TMe-717" 
#> [37] "TMe-2027" "TMe-3127" "TMe-3359" "TMe-2916" "TMe-1221" "TMe-1568"
#> [43] "TMe-1218" "TMe-500"  "TMe-306"  "TMe-3319" "TMe-1117" "TMe-3282"
#> [49] "TMe-937"  "TMe-3475" "TMe-501"  "TMe-3518" "TMe-3396" "TMe-3465"
#> [55] "TMe-2462" "TMe-2996" "TMe-3389" "TMe-2103" "TMe-3433" "TMe-3478"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-369"  "TMe-3466" "TMe-2056" "TMe-3605" "TMe-86"   "TMe-1107"
#>  [7] "TMe-47"   "TMe-774"  "TMe-160"  "TMe-674"  "TMe-404"  "TMe-3141"
#> [13] "TMe-3766" "TMe-1242" "TMe-196"  "TMe-2329" "TMe-289"  "TMe-1385"
#> [19] "TMe-2204" "TMe-1907" "TMe-3771" "TMe-176"  "TMe-1739" "TMe-3338"
#> [25] "TMe-478"  "TMe-509"  "TMe-3406" "TMe-3009" "TMe-2860" "TMe-2352"
#> [31] "TMe-171"  "TMe-3467" "TMe-2891" "TMe-477"  "TMe-3200" "TMe-400" 
#> [37] "TMe-539"  "TMe-1505" "TMe-2705" "TMe-3800" "TMe-3114"
#> 
#> $III
#>  [1] "TMe-878"  "TMe-3596" "TMe-785"  "TMe-2926" "TMe-3128" "TMe-1910"
#>  [7] "TMe-1897" "TMe-635"  "TMe-1863" "TMe-267"  "TMe-200"  "TMe-13"  
#> [13] "TMe-1443" "TMe-1804" "TMe-773"  "TMe-1889" "TMe-3274" "TMe-1790"
#> [19] "TMe-3383" "TMe-2823" "TMe-2394" "TMe-2811" "TMe-3565" "TMe-2481"
#> [25] "TMe-3352" "TMe-2502" "TMe-1939" "TMe-3556" "TMe-1797" "TMe-1868"
#> [31] "TMe-3216" "TMe-3234" "TMe-3445" "TMe-261"  "TMe-1838" "TMe-234" 
#> [37] "TMe-35"  
#> 
#> $IV
#>  [1] "TMe-2809" "TMe-609"  "TMe-1335" "TMe-241"  "TMe-596"  "TMe-1179"
#>  [7] "TMe-280"  "TMe-1975" "TMe-421"  "TMe-3417" "TMe-1525" "TMe-3386"
#> [13] "TMe-1988" "TMe-2788" "TMe-3040" "TMe-1095" "TMe-427"  "TMe-2833"
#> [19] "TMe-2956" "TMe-1434" "TMe-3025" "TMe-802"  "TMe-226"  "TMe-761" 
#> [25] "TMe-2567" "TMe-3273" "TMe-1336" "TMe-608"  "TMe-812"  "TMe-3316"
#> [31] "TMe-1078" "TMe-3581" "TMe-3576" "TMe-2043" "TMe-2958" "TMe-3758"
#> [37] "TMe-2924" "TMe-388"  "TMe-3055" "TMe-3442" "TMe-235"  "TMe-180" 
#> [43] "TMe-2960" "TMe-3054" "TMe-1129" "TMe-1891" "TMe-353"  "TMe-3701"
#> [49] "TMe-394"  "TMe-794"  "TMe-3752" "TMe-3573" "TMe-1903" "TMe-550" 
#> [55] "TMe-3730" "TMe-2059" "TMe-619"  "TMe-2318" "TMe-2928" "TMe-3017"
#> [61] "TMe-1923" "TMe-1016" "TMe-1376" "TMe-3253" "TMe-1419" "TMe-209" 
#> [67] "TMe-270"  "TMe-1231" "TMe-3594" "TMe-656"  "TMe-210"  "TMe-3072"
#> [73] "TMe-1597" "TMe-184"  "TMe-1297" "TMe-352"  "TMe-1020" "TMe-27"  
#> [79] "TMe-3527" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-399"  "TMe-3034" "TMe-2531" "TMe-277"  "TMe-1534" "TMe-700" 
#>  [7] "TMe-2121" "TMe-294"  "TMe-3628" "TMe-378"  "TMe-2853" "TMe-419" 
#> [13] "TMe-603"  "TMe-2933" "TMe-2050" "TMe-2907" "TMe-627"  "TMe-2296"
#> [19] "TMe-832"  "TMe-3411" "TMe-1311" "TMe-1399" "TMe-2790" "TMe-1011"
#> [25] "TMe-2589" "TMe-1307" "TMe-2003" "TMe-1294" "TMe-731"  "TMe-1158"
#> [31] "TMe-336"  "TMe-1099" "TMe-1455" "TMe-1257" "TMe-1300" "TMe-2435"
#> [37] "TMe-2953" "TMe-1730" "TMe-208"  "TMe-616"  "TMe-3736" "TMe-2290"
#> [43] "TMe-1762" "TMe-2688" "TMe-2761" "TMe-347"  "TMe-1290" "TMe-712" 
#> [49] "TMe-1458" "TMe-2355" "TMe-224"  "TMe-3408" "TMe-334"  "TMe-819" 
#> [55] "TMe-1372" "TMe-287"  "TMe-1391" "TMe-3659" "TMe-2750" "TMe-1269"
#> [61] "TMe-861"  "TMe-536"  "TMe-2060" "TMe-2820" "TMe-2413" "TMe-481" 
#> [67] "TMe-1501" "TMe-1488" "TMe-2119" "TMe-755"  "TMe-284"  "TMe-1446"
#> [73] "TMe-1098" "TMe-1924" "TMe-929"  "TMe-3332" "TMe-1366" "TMe-1381"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1428" "TMe-2963" "TMe-2035" "TMe-1232" "TMe-1403" "TMe-2983"
#>  [7] "TMe-1079" "TMe-3549" "TMe-1566" "TMe-1900" "TMe-3095" "TMe-705" 
#> [13] "TMe-1646" "TMe-1858" "TMe-836"  "TMe-3387" "TMe-690"  "TMe-893" 
#> [19] "TMe-1007" "TMe-2791" "TMe-1633" "TMe-1261" "TMe-725"  "TMe-1775"
#> [25] "TMe-2551" "TMe-1945" "TMe-1518" "TMe-1174" "TMe-1025" "TMe-1503"
#> [31] "TMe-1769" "TMe-2045" "TMe-856"  "TMe-1985" "TMe-1416" "TMe-1816"
#> [37] "TMe-3116"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_eccentricity_out, use.names = FALSE)) +
  labs(title = "eccentricity")


# Peripheral Sampling by Maximal Farness Centrality
sel_far_cent_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "farness.centrality")
sel_far_cent_out
#> $I
#>  [1] "TMe-2967" "TMe-3685" "TMe-3705" "TMe-3292" "TMe-3223" "TMe-3667"
#>  [7] "TMe-41"   "TMe-3475" "TMe-606"  "TMe-1425" "TMe-3282" "TMe-3319"
#> [13] "TMe-3392" "TMe-2965" "TMe-3437" "TMe-1564" "TMe-2943" "TMe-2955"
#> [19] "TMe-2993" "TMe-756"  "TMe-2069" "TMe-717"  "TMe-3415" "TMe-3396"
#> [25] "TMe-3389" "TMe-500"  "TMe-3065" "TMe-2996" "TMe-867"  "TMe-3641"
#> [31] "TMe-2779" "TMe-3115" "TMe-3163" "TMe-841"  "TMe-3151" "TMe-3433"
#> [37] "TMe-3249" "TMe-1360" "TMe-1190" "TMe-3330" "TMe-410"  "TMe-3266"
#> [43] "TMe-2027" "TMe-3102" "TMe-3398" "TMe-3264" "TMe-2604" "TMe-815" 
#> [49] "TMe-2785" "TMe-3478" "TMe-1218" "TMe-727"  "TMe-1333" "TMe-3323"
#> [55] "TMe-569"  "TMe-3471" "TMe-2985" "TMe-2984" "TMe-3281" "TMe-1466"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3466" "TMe-3338" "TMe-47"   "TMe-2860" "TMe-1919" "TMe-3605"
#>  [7] "TMe-1907" "TMe-86"   "TMe-1107" "TMe-2128" "TMe-369"  "TMe-2033"
#> [13] "TMe-2952" "TMe-2056" "TMe-3141" "TMe-3805" "TMe-2568" "TMe-2329"
#> [19] "TMe-1385" "TMe-1137" "TMe-3766" "TMe-3800" "TMe-674"  "TMe-3210"
#> [25] "TMe-3200" "TMe-2352" "TMe-2891" "TMe-3690" "TMe-3140" "TMe-509" 
#> [31] "TMe-1619" "TMe-1444" "TMe-478"  "TMe-404"  "TMe-2866" "TMe-171" 
#> [37] "TMe-289"  "TMe-2204" "TMe-160"  "TMe-477"  "TMe-3222"
#> 
#> $III
#>  [1] "TMe-3596" "TMe-2926" "TMe-1897" "TMe-3274" "TMe-3565" "TMe-267" 
#>  [7] "TMe-3804" "TMe-3299" "TMe-785"  "TMe-234"  "TMe-635"  "TMe-261" 
#> [13] "TMe-1868" "TMe-3445" "TMe-35"   "TMe-3663" "TMe-116"  "TMe-773" 
#> [19] "TMe-425"  "TMe-3631" "TMe-13"   "TMe-381"  "TMe-70"   "TMe-187" 
#> [25] "TMe-3032" "TMe-2823" "TMe-3383" "TMe-2897" "TMe-4"    "TMe-141" 
#> [31] "TMe-3143" "TMe-3043" "TMe-3715" "TMe-946"  "TMe-3048" "TMe-2756"
#> [37] "TMe-2811"
#> 
#> $IV
#>  [1] "TMe-2809" "TMe-241"  "TMe-1434" "TMe-550"  "TMe-812"  "TMe-761" 
#>  [7] "TMe-3417" "TMe-2924" "TMe-3442" "TMe-1988" "TMe-1525" "TMe-698" 
#> [13] "TMe-1020" "TMe-3730" "TMe-1078" "TMe-3054" "TMe-608"  "TMe-3273"
#> [19] "TMe-3573" "TMe-2043" "TMe-226"  "TMe-180"  "TMe-2956" "TMe-3275"
#> [25] "TMe-2567" "TMe-154"  "TMe-3545" "TMe-421"  "TMe-1229" "TMe-2971"
#> [31] "TMe-317"  "TMe-3025" "TMe-3707" "TMe-650"  "TMe-3167" "TMe-1162"
#> [37] "TMe-3189" "TMe-3255" "TMe-1129" "TMe-3297" "TMe-27"   "TMe-427" 
#> [43] "TMe-1419" "TMe-3256" "TMe-3055" "TMe-3089" "TMe-2958" "TMe-766" 
#> [49] "TMe-388"  "TMe-372"  "TMe-737"  "TMe-2979" "TMe-225"  "TMe-3701"
#> [55] "TMe-684"  "TMe-280"  "TMe-2552" "TMe-1376" "TMe-2240" "TMe-3225"
#> [61] "TMe-1328" "TMe-209"  "TMe-1479" "TMe-3053" "TMe-1377" "TMe-259" 
#> [67] "TMe-3594" "TMe-1348" "TMe-3527" "TMe-2833" "TMe-3040" "TMe-619" 
#> [73] "TMe-824"  "TMe-656"  "TMe-1891" "TMe-962"  "TMe-540"  "TMe-184" 
#> [79] "TMe-787"  "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-3034" "TMe-2121" "TMe-3628" "TMe-419"  "TMe-3736" "TMe-3185"
#>  [7] "TMe-536"  "TMe-2853" "TMe-378"  "TMe-712"  "TMe-2355" "TMe-2589"
#> [13] "TMe-2050" "TMe-1458" "TMe-2531" "TMe-616"  "TMe-723"  "TMe-3356"
#> [19] "TMe-1381" "TMe-1099" "TMe-1455" "TMe-3332" "TMe-399"  "TMe-700" 
#> [25] "TMe-782"  "TMe-2612" "TMe-603"  "TMe-277"  "TMe-861"  "TMe-2953"
#> [31] "TMe-284"  "TMe-1391" "TMe-2907" "TMe-2820" "TMe-2296" "TMe-3411"
#> [37] "TMe-2413" "TMe-3408" "TMe-2933" "TMe-2790" "TMe-832"  "TMe-2973"
#> [43] "TMe-2435" "TMe-3363" "TMe-1257" "TMe-1399" "TMe-1730" "TMe-3311"
#> [49] "TMe-2003" "TMe-344"  "TMe-2905" "TMe-2060" "TMe-189"  "TMe-1042"
#> [55] "TMe-2290" "TMe-294"  "TMe-2688" "TMe-1924" "TMe-48"   "TMe-969" 
#> [61] "TMe-1357" "TMe-1159" "TMe-2750" "TMe-730"  "TMe-208"  "TMe-532" 
#> [67] "TMe-3659" "TMe-2139" "TMe-1311" "TMe-2761" "TMe-755"  "TMe-1963"
#> [73] "TMe-158"  "TMe-481"  "TMe-2863" "TMe-287"  "TMe-2904" "TMe-2307"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-3095" "TMe-3549" "TMe-2983" "TMe-2035" "TMe-3116" "TMe-1416"
#>  [7] "TMe-693"  "TMe-2963" "TMe-1232" "TMe-1769" "TMe-1403" "TMe-2064"
#> [13] "TMe-1875" "TMe-1124" "TMe-1383" "TMe-1518" "TMe-1428" "TMe-1174"
#> [19] "TMe-2510" "TMe-1775" "TMe-725"  "TMe-1509" "TMe-1992" "TMe-465" 
#> [25] "TMe-1035" "TMe-1985" "TMe-2791" "TMe-836"  "TMe-690"  "TMe-1261"
#> [31] "TMe-1182" "TMe-3387" "TMe-1506" "TMe-751"  "TMe-2534" "TMe-361" 
#> [37] "TMe-1053"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_far_cent_out, use.names = FALSE)) +
  labs(title = "farness.centrality")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by space-Filling/coverage methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Space-Filling Sampling via the Kennard-Stone Algorithm
sel_ks_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "kennard.stone")
sel_ks_out
#> $I
#>  [1] "TMe-3268" "TMe-937"  "TMe-3623" "TMe-1621" "TMe-3481" "TMe-3719"
#>  [7] "TMe-501"  "TMe-2984" "TMe-3163" "TMe-3208" "TMe-469"  "TMe-815" 
#> [13] "TMe-1568" "TMe-3330" "TMe-3521" "TMe-438"  "TMe-1347" "TMe-3145"
#> [19] "TMe-3219" "TMe-2462" "TMe-3132" "TMe-3292" "TMe-3392" "TMe-3518"
#> [25] "TMe-33"   "TMe-446"  "TMe-736"  "TMe-778"  "TMe-867"  "TMe-1096"
#> [31] "TMe-1448" "TMe-1823" "TMe-2031" "TMe-2255" "TMe-3027" "TMe-3641"
#> [37] "TMe-1532" "TMe-1589" "TMe-2909" "TMe-2955" "TMe-3164" "TMe-3351"
#> [43] "TMe-3398" "TMe-3419" "TMe-44"   "TMe-290"  "TMe-489"  "TMe-940" 
#> [49] "TMe-1191" "TMe-1341" "TMe-1466" "TMe-1906" "TMe-1930" "TMe-2191"
#> [55] "TMe-2785" "TMe-3030" "TMe-3127" "TMe-3272" "TMe-2912" "TMe-3001"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3571" "TMe-589"  "TMe-890"  "TMe-660"  "TMe-1409" "TMe-706" 
#>  [7] "TMe-539"  "TMe-1831" "TMe-3447" "TMe-3467" "TMe-40"   "TMe-1619"
#> [13] "TMe-796"  "TMe-2200" "TMe-2851" "TMe-2978" "TMe-3009" "TMe-3222"
#> [19] "TMe-774"  "TMe-2"    "TMe-85"   "TMe-1051" "TMe-1754" "TMe-1864"
#> [25] "TMe-2304" "TMe-2797" "TMe-3114" "TMe-3557" "TMe-196"  "TMe-1459"
#> [31] "TMe-1616" "TMe-160"  "TMe-377"  "TMe-477"  "TMe-1005" "TMe-1137"
#> [37] "TMe-1250" "TMe-1732" "TMe-2123" "TMe-2862" "TMe-3052"
#> 
#> $III
#>  [1] "TMe-3715" "TMe-1421" "TMe-3287" "TMe-1593" "TMe-14"   "TMe-434" 
#>  [7] "TMe-141"  "TMe-70"   "TMe-2926" "TMe-3008" "TMe-773"  "TMe-2270"
#> [13] "TMe-203"  "TMe-381"  "TMe-2086" "TMe-2502" "TMe-1230" "TMe-425" 
#> [19] "TMe-1819" "TMe-3148" "TMe-6"    "TMe-261"  "TMe-785"  "TMe-1176"
#> [25] "TMe-1200" "TMe-1790" "TMe-234"  "TMe-1207" "TMe-3317" "TMe-3336"
#> [31] "TMe-206"  "TMe-1838" "TMe-2751" "TMe-3362" "TMe-3663" "TMe-3679"
#> [37] "TMe-3731"
#> 
#> $IV
#>  [1] "TMe-3340" "TMe-90"   "TMe-1776" "TMe-519"  "TMe-27"   "TMe-12"  
#>  [7] "TMe-2458" "TMe-1526" "TMe-3144" "TMe-3253" "TMe-107"  "TMe-526" 
#> [13] "TMe-2954" "TMe-1419" "TMe-3276" "TMe-57"   "TMe-204"  "TMe-1095"
#> [19] "TMe-1155" "TMe-2764" "TMe-3191" "TMe-81"   "TMe-887"  "TMe-3412"
#> [25] "TMe-3639" "TMe-15"   "TMe-3484" "TMe-154"  "TMe-396"  "TMe-760" 
#> [31] "TMe-1246" "TMe-1814" "TMe-2039" "TMe-3072" "TMe-3435" "TMe-30"  
#> [37] "TMe-153"  "TMe-207"  "TMe-280"  "TMe-1651" "TMe-2989" "TMe-3290"
#> [43] "TMe-3494" "TMe-3558" "TMe-3214" "TMe-534"  "TMe-812"  "TMe-840" 
#> [49] "TMe-1351" "TMe-1525" "TMe-1700" "TMe-1923" "TMe-3259" "TMe-3357"
#> [55] "TMe-7"    "TMe-67"   "TMe-1150" "TMe-3019" "TMe-3167" "TMe-3701"
#> [61] "TMe-2960" "TMe-266"  "TMe-403"  "TMe-525"  "TMe-641"  "TMe-824" 
#> [67] "TMe-919"  "TMe-1166" "TMe-1348" "TMe-1364" "TMe-2809" "TMe-2971"
#> [73] "TMe-3068" "TMe-3103" "TMe-3161" "TMe-3257" "TMe-3273" "TMe-25"  
#> [79] "TMe-87"   "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-1551" "TMe-373"  "TMe-600"  "TMe-1345" "TMe-1733" "TMe-1411"
#>  [7] "TMe-487"  "TMe-1534" "TMe-886"  "TMe-1160" "TMe-1418" "TMe-2213"
#> [13] "TMe-167"  "TMe-217"  "TMe-655"  "TMe-850"  "TMe-399"  "TMe-997" 
#> [19] "TMe-1500" "TMe-419"  "TMe-456"  "TMe-1268" "TMe-1373" "TMe-1880"
#> [25] "TMe-2121" "TMe-1322" "TMe-347"  "TMe-378"  "TMe-2530" "TMe-168" 
#> [31] "TMe-712"  "TMe-1788" "TMe-3311" "TMe-423"  "TMe-474"  "TMe-603" 
#> [37] "TMe-1037" "TMe-1963" "TMe-476"  "TMe-532"  "TMe-945"  "TMe-1381"
#> [43] "TMe-2138" "TMe-771"  "TMe-927"  "TMe-1009" "TMe-1762" "TMe-3332"
#> [49] "TMe-892"  "TMe-574"  "TMe-629"  "TMe-777"  "TMe-1257" "TMe-1440"
#> [55] "TMe-1862" "TMe-1901" "TMe-2435" "TMe-256"  "TMe-627"  "TMe-1158"
#> [61] "TMe-1283" "TMe-1311" "TMe-307"  "TMe-414"  "TMe-679"  "TMe-797" 
#> [67] "TMe-1204" "TMe-1834" "TMe-1871" "TMe-2060" "TMe-2589" "TMe-2821"
#> [73] "TMe-2863" "TMe-193"  "TMe-390"  "TMe-464"  "TMe-584"  "TMe-861" 
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1509" "TMe-1985" "TMe-1481" "TMe-693"  "TMe-1233" "TMe-696" 
#>  [7] "TMe-2791" "TMe-222"  "TMe-3549" "TMe-631"  "TMe-1483" "TMe-531" 
#> [13] "TMe-985"  "TMe-1261" "TMe-1413" "TMe-1428" "TMe-662"  "TMe-1178"
#> [19] "TMe-1678" "TMe-1744" "TMe-315"  "TMe-361"  "TMe-1518" "TMe-1775"
#> [25] "TMe-2196" "TMe-213"  "TMe-728"  "TMe-1053" "TMe-1506" "TMe-1723"
#> [31] "TMe-3095" "TMe-321"  "TMe-725"  "TMe-845"  "TMe-1362" "TMe-1445"
#> [37] "TMe-2045"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_ks_out, use.names = FALSE)) +
  labs(title = "kennard.stone")


# Space-Filling Sampling via the DUPLEX Algorithm
sel_duplex_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "duplex")
sel_duplex_out
#> $I
#>  [1] "TMe-3268" "TMe-937"  "TMe-3623" "TMe-1621" "TMe-3481" "TMe-3719"
#>  [7] "TMe-501"  "TMe-2984" "TMe-3163" "TMe-3208" "TMe-469"  "TMe-1568"
#> [13] "TMe-3330" "TMe-438"  "TMe-1347" "TMe-3145" "TMe-2462" "TMe-2909"
#> [19] "TMe-3132" "TMe-3292" "TMe-3518" "TMe-33"   "TMe-446"  "TMe-736" 
#> [25] "TMe-778"  "TMe-867"  "TMe-1096" "TMe-1448" "TMe-1823" "TMe-2388"
#> [31] "TMe-3027" "TMe-1532" "TMe-1589" "TMe-2955" "TMe-3398" "TMe-44"  
#> [37] "TMe-290"  "TMe-1140" "TMe-1191" "TMe-1830" "TMe-1906" "TMe-1930"
#> [43] "TMe-2785" "TMe-3127" "TMe-3112" "TMe-3142" "TMe-3252" "TMe-3281"
#> [49] "TMe-3325" "TMe-3429" "TMe-566"  "TMe-642"  "TMe-1156" "TMe-1869"
#> [55] "TMe-1960" "TMe-1964" "TMe-2004" "TMe-2152" "TMe-2604" "TMe-28"  
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3571" "TMe-589"  "TMe-890"  "TMe-660"  "TMe-1409" "TMe-706" 
#>  [7] "TMe-539"  "TMe-2715" "TMe-3447" "TMe-40"   "TMe-1619" "TMe-796" 
#> [13] "TMe-1137" "TMe-2851" "TMe-85"   "TMe-339"  "TMe-377"  "TMe-1051"
#> [19] "TMe-1250" "TMe-1732" "TMe-2304" "TMe-2797" "TMe-1459" "TMe-3443"
#> [25] "TMe-2"    "TMe-369"  "TMe-1005" "TMe-1754" "TMe-2054" "TMe-2862"
#> [31] "TMe-53"   "TMe-821"  "TMe-1251" "TMe-1310" "TMe-1668" "TMe-2843"
#> [37] "TMe-3599" "TMe-3690" "TMe-289"  "TMe-365"  "TMe-1279"
#> 
#> $III
#>  [1] "TMe-3715" "TMe-1421" "TMe-3287" "TMe-1593" "TMe-14"   "TMe-2926"
#>  [7] "TMe-141"  "TMe-2203" "TMe-1910" "TMe-6"    "TMe-206"  "TMe-3362"
#> [13] "TMe-70"   "TMe-3143" "TMe-1176" "TMe-1868" "TMe-2161" "TMe-3147"
#> [19] "TMe-2167" "TMe-2751" "TMe-3796" "TMe-2394" "TMe-2977" "TMe-3569"
#> [25] "TMe-3700" "TMe-425"  "TMe-773"  "TMe-1797" "TMe-2481" "TMe-3346"
#> [31] "TMe-3679" "TMe-3731" "TMe-3383" "TMe-420"  "TMe-1817" "TMe-1819"
#> [37] "TMe-2205"
#> 
#> $IV
#>  [1] "TMe-3340" "TMe-90"   "TMe-1776" "TMe-519"  "TMe-27"   "TMe-12"  
#>  [7] "TMe-2458" "TMe-1526" "TMe-3144" "TMe-3253" "TMe-526"  "TMe-2954"
#> [13] "TMe-1419" "TMe-1765" "TMe-3276" "TMe-204"  "TMe-834"  "TMe-1155"
#> [19] "TMe-3257" "TMe-887"  "TMe-3639" "TMe-3484" "TMe-396"  "TMe-760" 
#> [25] "TMe-1246" "TMe-3072" "TMe-3435" "TMe-30"   "TMe-207"  "TMe-1651"
#> [31] "TMe-2989" "TMe-3019" "TMe-3494" "TMe-3558" "TMe-15"   "TMe-840" 
#> [37] "TMe-1123" "TMe-1351" "TMe-1806" "TMe-1923" "TMe-3259" "TMe-3779"
#> [43] "TMe-67"   "TMe-87"   "TMe-421"  "TMe-634"  "TMe-1336" "TMe-2375"
#> [49] "TMe-919"  "TMe-1368" "TMe-1948" "TMe-2807" "TMe-2971" "TMe-3073"
#> [55] "TMe-3440" "TMe-180"  "TMe-708"  "TMe-761"  "TMe-1016" "TMe-1328"
#> [61] "TMe-1456" "TMe-1479" "TMe-3040" "TMe-3109" "TMe-3225" "TMe-3660"
#> [67] "TMe-191"  "TMe-352"  "TMe-675"  "TMe-1118" "TMe-1151" "TMe-1485"
#> [73] "TMe-1867" "TMe-2043" "TMe-2981" "TMe-3255" "TMe-54"   "TMe-72"  
#> [79] "TMe-82"   "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-1551" "TMe-373"  "TMe-600"  "TMe-1345" "TMe-1733" "TMe-1411"
#>  [7] "TMe-487"  "TMe-1534" "TMe-1160" "TMe-1418" "TMe-2213" "TMe-217" 
#> [13] "TMe-655"  "TMe-850"  "TMe-399"  "TMe-2015" "TMe-419"  "TMe-1268"
#> [19] "TMe-1373" "TMe-1880" "TMe-748"  "TMe-1322" "TMe-347"  "TMe-2530"
#> [25] "TMe-168"  "TMe-370"  "TMe-712"  "TMe-3311" "TMe-474"  "TMe-603" 
#> [31] "TMe-667"  "TMe-1963" "TMe-771"  "TMe-532"  "TMe-574"  "TMe-945" 
#> [37] "TMe-1381" "TMe-1762" "TMe-2953" "TMe-3332" "TMe-629"  "TMe-777" 
#> [43] "TMe-1440" "TMe-1862" "TMe-1871" "TMe-2435" "TMe-211"  "TMe-481" 
#> [49] "TMe-1158" "TMe-414"  "TMe-679"  "TMe-797"  "TMe-861"  "TMe-1204"
#> [55] "TMe-1257" "TMe-1834" "TMe-1934" "TMe-2590" "TMe-193"  "TMe-901" 
#> [61] "TMe-1042" "TMe-3408" "TMe-892"  "TMe-464"  "TMe-929"  "TMe-1357"
#> [67] "TMe-1427" "TMe-1873" "TMe-2151" "TMe-2589" "TMe-430"  "TMe-816" 
#> [73] "TMe-1788" "TMe-3356" "TMe-682"  "TMe-1859" "TMe-247"  "TMe-647" 
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1509" "TMe-1985" "TMe-1481" "TMe-693"  "TMe-1233" "TMe-696" 
#>  [7] "TMe-2791" "TMe-3549" "TMe-631"  "TMe-1483" "TMe-662"  "TMe-531" 
#> [13] "TMe-985"  "TMe-1261" "TMe-1413" "TMe-1428" "TMe-1678" "TMe-1744"
#> [19] "TMe-315"  "TMe-1775" "TMe-188"  "TMe-213"  "TMe-1506" "TMe-1723"
#> [25] "TMe-2535" "TMe-3095" "TMe-321"  "TMe-725"  "TMe-1362" "TMe-1445"
#> [31] "TMe-1992" "TMe-2045" "TMe-852"  "TMe-963"  "TMe-1079" "TMe-1319"
#> [37] "TMe-1477"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_duplex_out, use.names = FALSE)) +
  labs(title = "duplex")


# Space-Filling Sampling via the Honigs Algorithm
sel_honigs_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "honigs")
sel_honigs_out
#> $I
#>  [1] "TMe-3667" "TMe-606"  "TMe-2965" "TMe-1423" "TMe-2967" "TMe-3337"
#>  [7] "TMe-3223" "TMe-2513" "TMe-3163" "TMe-3151" "TMe-3478" "TMe-2993"
#> [13] "TMe-2934" "TMe-1786" "TMe-3481" "TMe-1589" "TMe-715"  "TMe-264" 
#> [19] "TMe-1425" "TMe-867"  "TMe-2975" "TMe-3026" "TMe-299"  "TMe-3437"
#> [25] "TMe-579"  "TMe-3471" "TMe-3175" "TMe-132"  "TMe-2462" "TMe-3292"
#> [31] "TMe-3460" "TMe-3463" "TMe-3030" "TMe-41"   "TMe-778"  "TMe-2027"
#> [37] "TMe-3425" "TMe-3353" "TMe-3142" "TMe-1960" "TMe-2010" "TMe-1466"
#> [43] "TMe-3396" "TMe-2031" "TMe-3102" "TMe-2955" "TMe-3415" "TMe-2985"
#> [49] "TMe-3169" "TMe-3601" "TMe-717"  "TMe-3685" "TMe-610"  "TMe-815" 
#> [55] "TMe-1564" "TMe-937"  "TMe-2913" "TMe-1226" "TMe-2906" "TMe-1940"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3466" "TMe-369"  "TMe-3605" "TMe-1919" "TMe-1907" "TMe-3805"
#>  [7] "TMe-47"   "TMe-171"  "TMe-3338" "TMe-289"  "TMe-2329" "TMe-3690"
#> [13] "TMe-3366" "TMe-3141" "TMe-3547" "TMe-431"  "TMe-3150" "TMe-160" 
#> [19] "TMe-3800" "TMe-3209" "TMe-1051" "TMe-2204" "TMe-3530" "TMe-2891"
#> [25] "TMe-2952" "TMe-60"   "TMe-2033" "TMe-3766" "TMe-3222" "TMe-3286"
#> [31] "TMe-509"  "TMe-1668" "TMe-2866" "TMe-902"  "TMe-1461" "TMe-2304"
#> [37] "TMe-3771" "TMe-2860" "TMe-1310" "TMe-2056" "TMe-1137"
#> 
#> $III
#>  [1] "TMe-3596" "TMe-878"  "TMe-234"  "TMe-425"  "TMe-1897" "TMe-946" 
#>  [7] "TMe-35"   "TMe-3274" "TMe-13"   "TMe-6"    "TMe-1993" "TMe-3299"
#> [13] "TMe-1200" "TMe-3679" "TMe-3715" "TMe-2374" "TMe-3700" "TMe-2086"
#> [19] "TMe-785"  "TMe-3383" "TMe-267"  "TMe-3575" "TMe-2756" "TMe-3048"
#> [25] "TMe-2733" "TMe-3804" "TMe-2823" "TMe-1868" "TMe-2592" "TMe-3608"
#> [31] "TMe-2926" "TMe-3397" "TMe-2270" "TMe-3565" "TMe-2977" "TMe-3407"
#> [37] "TMe-1176"
#> 
#> $IV
#>  [1] "TMe-609"  "TMe-2809" "TMe-761"  "TMe-3054" "TMe-550"  "TMe-812" 
#>  [7] "TMe-2924" "TMe-225"  "TMe-2971" "TMe-3025" "TMe-3542" "TMe-3273"
#> [13] "TMe-3422" "TMe-656"  "TMe-1020" "TMe-3428" "TMe-3527" "TMe-3275"
#> [19] "TMe-2240" "TMe-3573" "TMe-3707" "TMe-3730" "TMe-317"  "TMe-226" 
#> [25] "TMe-3593" "TMe-2567" "TMe-733"  "TMe-372"  "TMe-241"  "TMe-1511"
#> [31] "TMe-3055" "TMe-1988" "TMe-919"  "TMe-3248" "TMe-2039" "TMe-1397"
#> [37] "TMe-3390" "TMe-897"  "TMe-3701" "TMe-3255" "TMe-1575" "TMe-2043"
#> [43] "TMe-698"  "TMe-1083" "TMe-427"  "TMe-388"  "TMe-3619" "TMe-3297"
#> [49] "TMe-1700" "TMe-2981" "TMe-650"  "TMe-1376" "TMe-57"   "TMe-3386"
#> [55] "TMe-2776" "TMe-737"  "TMe-3089" "TMe-1806" "TMe-180"  "TMe-3167"
#> [61] "TMe-66"   "TMe-3040" "TMe-556"  "TMe-526"  "TMe-3004" "TMe-5"   
#> [67] "TMe-3417" "TMe-1434" "TMe-18"   "TMe-1348" "TMe-1123" "TMe-2025"
#> [73] "TMe-1297" "TMe-540"  "TMe-72"   "TMe-3340" "TMe-1148" "TMe-3103"
#> [79] "TMe-3545" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-3034" "TMe-399"  "TMe-712"  "TMe-3185" "TMe-1381" "TMe-3332"
#>  [7] "TMe-277"  "TMe-2853" "TMe-2050" "TMe-832"  "TMe-723"  "TMe-700" 
#> [13] "TMe-2933" "TMe-798"  "TMe-1366" "TMe-3736" "TMe-870"  "TMe-759" 
#> [19] "TMe-487"  "TMe-895"  "TMe-336"  "TMe-1290" "TMe-3628" "TMe-1401"
#> [25] "TMe-208"  "TMe-2413" "TMe-2790" "TMe-1042" "TMe-2905" "TMe-2296"
#> [31] "TMe-2688" "TMe-1026" "TMe-782"  "TMe-2139" "TMe-2973" "TMe-391" 
#> [37] "TMe-536"  "TMe-3356" "TMe-1004" "TMe-588"  "TMe-2041" "TMe-418" 
#> [43] "TMe-1482" "TMe-1924" "TMe-2121" "TMe-2589" "TMe-943"  "TMe-861" 
#> [49] "TMe-819"  "TMe-2953" "TMe-603"  "TMe-334"  "TMe-1257" "TMe-2750"
#> [55] "TMe-2907" "TMe-1099" "TMe-2290" "TMe-863"  "TMe-901"  "TMe-99"  
#> [61] "TMe-1269" "TMe-1963" "TMe-547"  "TMe-2531" "TMe-1293" "TMe-1196"
#> [67] "TMe-378"  "TMe-1283" "TMe-1307" "TMe-412"  "TMe-397"  "TMe-616" 
#> [73] "TMe-1158" "TMe-1311" "TMe-1501" "TMe-969"  "TMe-532"  "TMe-731" 
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-2963" "TMe-1428" "TMe-1875" "TMe-1035" "TMe-1232" "TMe-361" 
#>  [7] "TMe-3095" "TMe-2035" "TMe-1769" "TMe-693"  "TMe-315"  "TMe-1124"
#> [13] "TMe-751"  "TMe-1796" "TMe-531"  "TMe-2791" "TMe-1174" "TMe-598" 
#> [19] "TMe-442"  "TMe-1809" "TMe-138"  "TMe-1816" "TMe-1404" "TMe-1858"
#> [25] "TMe-1646" "TMe-1302" "TMe-1592" "TMe-2534" "TMe-1756" "TMe-1383"
#> [31] "TMe-1518" "TMe-1025" "TMe-236"  "TMe-1416" "TMe-2398" "TMe-1566"
#> [37] "TMe-1362"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_honigs_out, use.names = FALSE)) +
  labs(title = "honigs")


# Space-Filling Sampling via Farthest-Point (Max-Min) Algorithm
sel_far_pt_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "farthest.point")
sel_far_pt_out
#> $I
#>  [1] "TMe-2216" "TMe-3705" "TMe-2965" "TMe-756"  "TMe-264"  "TMe-3437"
#>  [7] "TMe-1425" "TMe-3266" "TMe-2912" "TMe-1589" "TMe-3415" "TMe-778" 
#> [13] "TMe-2934" "TMe-300"  "TMe-3514" "TMe-1786" "TMe-2943" "TMe-1922"
#> [19] "TMe-3475" "TMe-566"  "TMe-3462" "TMe-3163" "TMe-2993" "TMe-715" 
#> [25] "TMe-410"  "TMe-3478" "TMe-3292" "TMe-1958" "TMe-3625" "TMe-2955"
#> [31] "TMe-3208" "TMe-2027" "TMe-2939" "TMe-3223" "TMe-41"   "TMe-2975"
#> [37] "TMe-2937" "TMe-2462" "TMe-3457" "TMe-3553" "TMe-3389" "TMe-1466"
#> [43] "TMe-2779" "TMe-2916" "TMe-3249" "TMe-3499" "TMe-2785" "TMe-3026"
#> [49] "TMe-3127" "TMe-717"  "TMe-2010" "TMe-867"  "TMe-3351" "TMe-3515"
#> [55] "TMe-3051" "TMe-2372" "TMe-727"  "TMe-3087" "TMe-2945" "TMe-3169"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2352" "TMe-369"  "TMe-3605" "TMe-2056" "TMe-3466" "TMe-1459"
#>  [7] "TMe-2860" "TMe-3690" "TMe-47"   "TMe-3447" "TMe-3210" "TMe-3338"
#> [13] "TMe-2033" "TMe-509"  "TMe-1619" "TMe-768"  "TMe-1919" "TMe-3455"
#> [19] "TMe-3800" "TMe-3406" "TMe-3200" "TMe-1385" "TMe-902"  "TMe-289" 
#> [25] "TMe-3052" "TMe-2952" "TMe-3009" "TMe-648"  "TMe-1907" "TMe-1754"
#> [31] "TMe-2128" "TMe-2304" "TMe-1668" "TMe-2000" "TMe-2797" "TMe-339" 
#> [37] "TMe-925"  "TMe-2"    "TMe-3547" "TMe-1241" "TMe-3141"
#> 
#> $III
#>  [1] "TMe-2205" "TMe-3565" "TMe-3596" "TMe-635"  "TMe-2977" "TMe-1897"
#>  [7] "TMe-3274" "TMe-2203" "TMe-4"    "TMe-3299" "TMe-3721" "TMe-35"  
#> [13] "TMe-234"  "TMe-123"  "TMe-3445" "TMe-3715" "TMe-3397" "TMe-785" 
#> [19] "TMe-2901" "TMe-2897" "TMe-3048" "TMe-3230" "TMe-3796" "TMe-2481"
#> [25] "TMe-3143" "TMe-3383" "TMe-200"  "TMe-1939" "TMe-141"  "TMe-261" 
#> [31] "TMe-3551" "TMe-425"  "TMe-2592" "TMe-946"  "TMe-2161" "TMe-267" 
#> [37] "TMe-830" 
#> 
#> $IV
#>  [1] "TMe-3417" "TMe-421"  "TMe-3273" "TMe-698"  "TMe-1700" "TMe-2809"
#>  [7] "TMe-3701" "TMe-2240" "TMe-3054" "TMe-1434" "TMe-2775" "TMe-619" 
#> [13] "TMe-226"  "TMe-3573" "TMe-761"  "TMe-3730" "TMe-824"  "TMe-766" 
#> [19] "TMe-1350" "TMe-3760" "TMe-2839" "TMe-3576" "TMe-383"  "TMe-3297"
#> [25] "TMe-427"  "TMe-184"  "TMe-1278" "TMe-3428" "TMe-967"  "TMe-737" 
#> [31] "TMe-3422" "TMe-2956" "TMe-3055" "TMe-1806" "TMe-1867" "TMe-5"   
#> [37] "TMe-259"  "TMe-3312" "TMe-812"  "TMe-2043" "TMe-3189" "TMe-3780"
#> [43] "TMe-2025" "TMe-526"  "TMe-3225" "TMe-1148" "TMe-3283" "TMe-204" 
#> [49] "TMe-3276" "TMe-3025" "TMe-1368" "TMe-3256" "TMe-787"  "TMe-3390"
#> [55] "TMe-3440" "TMe-3017" "TMe-2979" "TMe-1597" "TMe-3214" "TMe-2776"
#> [61] "TMe-594"  "TMe-2552" "TMe-3089" "TMe-596"  "TMe-1814" "TMe-1525"
#> [67] "TMe-190"  "TMe-1325" "TMe-3040" "TMe-25"   "TMe-180"  "TMe-1348"
#> [73] "TMe-270"  "TMe-225"  "TMe-1511" "TMe-72"   "TMe-3255" "TMe-1988"
#> [79] "TMe-2981" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-2139" "TMe-700"  "TMe-2531" "TMe-3311" "TMe-3736" "TMe-2853"
#>  [7] "TMe-2863" "TMe-2121" "TMe-194"  "TMe-1500" "TMe-481"  "TMe-655" 
#> [13] "TMe-603"  "TMe-284"  "TMe-336"  "TMe-3185" "TMe-616"  "TMe-1257"
#> [19] "TMe-294"  "TMe-2589" "TMe-3034" "TMe-362"  "TMe-816"  "TMe-582" 
#> [25] "TMe-819"  "TMe-2907" "TMe-3628" "TMe-287"  "TMe-1381" "TMe-1311"
#> [31] "TMe-3411" "TMe-305"  "TMe-378"  "TMe-277"  "TMe-2355" "TMe-1643"
#> [37] "TMe-588"  "TMe-208"  "TMe-1301" "TMe-158"  "TMe-1026" "TMe-536" 
#> [43] "TMe-3329" "TMe-1802" "TMe-2060" "TMe-1196" "TMe-2973" "TMe-3408"
#> [49] "TMe-419"  "TMe-2905" "TMe-1730" "TMe-769"  "TMe-3659" "TMe-712" 
#> [55] "TMe-2790" "TMe-1204" "TMe-48"   "TMe-2959" "TMe-1401" "TMe-1188"
#> [61] "TMe-1372" "TMe-2413" "TMe-782"  "TMe-1042" "TMe-2688" "TMe-977" 
#> [67] "TMe-2192" "TMe-3356" "TMe-1963" "TMe-798"  "TMe-861"  "TMe-1366"
#> [73] "TMe-1009" "TMe-3363" "TMe-877"  "TMe-688"  "TMe-997"  "TMe-2425"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1403" "TMe-2035" "TMe-580"  "TMe-693"  "TMe-2534" "TMe-751" 
#>  [7] "TMe-1232" "TMe-1985" "TMe-2791" "TMe-1796" "TMe-2045" "TMe-1079"
#> [13] "TMe-1416" "TMe-310"  "TMe-548"  "TMe-1775" "TMe-1174" "TMe-1646"
#> [19] "TMe-1995" "TMe-690"  "TMe-1264" "TMe-2026" "TMe-361"  "TMe-130" 
#> [25] "TMe-598"  "TMe-1518" "TMe-1900" "TMe-1076" "TMe-1769" "TMe-1302"
#> [31] "TMe-1025" "TMe-1035" "TMe-696"  "TMe-2510" "TMe-836"  "TMe-3549"
#> [37] "TMe-2067"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_far_pt_out, use.names = FALSE)) +
  labs(title = "farthest.point")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by density based methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Density-Based Sampling by Minimal Nearest-Neighbour Distance
sel_nn_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "nearest.neighbour")
sel_nn_out
#> $I
#>  [1] "TMe-1696" "TMe-1717" "TMe-2976" "TMe-20"   "TMe-2913" "TMe-2916"
#>  [7] "TMe-839"  "TMe-1086" "TMe-1140" "TMe-2084" "TMe-2103" "TMe-2152"
#> [13] "TMe-1423" "TMe-2131" "TMe-438"  "TMe-3424" "TMe-1448" "TMe-2191"
#> [19] "TMe-1221" "TMe-3170" "TMe-3548" "TMe-486"  "TMe-2251" "TMe-990" 
#> [25] "TMe-1823" "TMe-3325" "TMe-3729" "TMe-583"  "TMe-1981" "TMe-3027"
#> [31] "TMe-1117" "TMe-910"  "TMe-3268" "TMe-3235" "TMe-1955" "TMe-1869"
#> [37] "TMe-2206" "TMe-3149" "TMe-3208" "TMe-2066" "TMe-2383" "TMe-501" 
#> [43] "TMe-2255" "TMe-248"  "TMe-686"  "TMe-1344" "TMe-3202" "TMe-2388"
#> [49] "TMe-3115" "TMe-3065" "TMe-290"  "TMe-2914" "TMe-865"  "TMe-1091"
#> [55] "TMe-1224" "TMe-579"  "TMe-1930" "TMe-933"  "TMe-489"  "TMe-2949"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-455"  "TMe-890"  "TMe-3066" "TMe-3571" "TMe-528"  "TMe-930" 
#>  [7] "TMe-1505" "TMe-2814" "TMe-2127" "TMe-664"  "TMe-1251" "TMe-1833"
#> [13] "TMe-3203" "TMe-1409" "TMe-1766" "TMe-2021" "TMe-2412" "TMe-365" 
#> [19] "TMe-2611" "TMe-3533" "TMe-377"  "TMe-2935" "TMe-433"  "TMe-2705"
#> [25] "TMe-2765" "TMe-3187" "TMe-2951" "TMe-1969" "TMe-3642" "TMe-3531"
#> [31] "TMe-3239" "TMe-3020" "TMe-2242" "TMe-3540" "TMe-1795" "TMe-60"  
#> [37] "TMe-1323" "TMe-3495" "TMe-3617" "TMe-39"   "TMe-589" 
#> 
#> $III
#>  [1] "TMe-2161" "TMe-10"   "TMe-853"  "TMe-1738" "TMe-304"  "TMe-1836"
#>  [7] "TMe-3631" "TMe-3721" "TMe-64"   "TMe-3346" "TMe-1797" "TMe-1910"
#> [13] "TMe-3118" "TMe-617"  "TMe-3113" "TMe-3147" "TMe-3336" "TMe-1790"
#> [19] "TMe-1804" "TMe-3005" "TMe-3335" "TMe-161"  "TMe-3731" "TMe-1138"
#> [25] "TMe-2748" "TMe-3324" "TMe-1207" "TMe-1675" "TMe-3069" "TMe-3544"
#> [31] "TMe-3572" "TMe-237"  "TMe-1965" "TMe-3043" "TMe-3148" "TMe-2167"
#> [37] "TMe-2154"
#> 
#> $IV
#>  [1] "TMe-2247" "TMe-1139" "TMe-3231" "TMe-3660" "TMe-1118" "TMe-1402"
#>  [7] "TMe-107"  "TMe-72"   "TMe-885"  "TMe-2410" "TMe-57"   "TMe-204" 
#> [13] "TMe-150"  "TMe-1364" "TMe-24"   "TMe-3044" "TMe-3243" "TMe-3779"
#> [19] "TMe-1726" "TMe-1948" "TMe-3357" "TMe-7"    "TMe-380"  "TMe-218" 
#> [25] "TMe-396"  "TMe-3639" "TMe-3494" "TMe-103"  "TMe-3233" "TMe-3240"
#> [31] "TMe-428"  "TMe-1765" "TMe-1374" "TMe-1470" "TMe-2210" "TMe-887" 
#> [37] "TMe-581"  "TMe-2775" "TMe-3073" "TMe-67"   "TMe-1377" "TMe-1419"
#> [43] "TMe-875"  "TMe-1652" "TMe-3232" "TMe-1106" "TMe-1123" "TMe-2039"
#> [49] "TMe-3562" "TMe-3412" "TMe-3757" "TMe-12"   "TMe-3238" "TMe-3084"
#> [55] "TMe-3780" "TMe-1166" "TMe-840"  "TMe-3378" "TMe-1211" "TMe-3161"
#> [61] "TMe-1167" "TMe-824"  "TMe-834"  "TMe-54"   "TMe-3045" "TMe-2340"
#> [67] "TMe-2758" "TMe-3454" "TMe-3198" "TMe-3204" "TMe-641"  "TMe-3781"
#> [73] "TMe-3615" "TMe-526"  "TMe-3435" "TMe-467"  "TMe-3168" "TMe-108" 
#> [79] "TMe-525"  "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-424"  "TMe-430"  "TMe-1559" "TMe-977"  "TMe-1234" "TMe-1572"
#>  [7] "TMe-1879" "TMe-2009" "TMe-2441" "TMe-2361" "TMe-212"  "TMe-464" 
#> [13] "TMe-142"  "TMe-669"  "TMe-1290" "TMe-474"  "TMe-1414" "TMe-618" 
#> [19] "TMe-889"  "TMe-1160" "TMe-1862" "TMe-137"  "TMe-360"  "TMe-1521"
#> [25] "TMe-565"  "TMe-2845" "TMe-168"  "TMe-1188" "TMe-754"  "TMe-330" 
#> [31] "TMe-1367" "TMe-1901" "TMe-1294" "TMe-2439" "TMe-1871" "TMe-1541"
#> [37] "TMe-2213" "TMe-2319" "TMe-1979" "TMe-343"  "TMe-382"  "TMe-939" 
#> [43] "TMe-363"  "TMe-473"  "TMe-886"  "TMe-1011" "TMe-822"  "TMe-2124"
#> [49] "TMe-685"  "TMe-1308" "TMe-312"  "TMe-1733" "TMe-1778" "TMe-307" 
#> [55] "TMe-348"  "TMe-647"  "TMe-1629" "TMe-1312" "TMe-1382" "TMe-456" 
#> [61] "TMe-835"  "TMe-510"  "TMe-870"  "TMe-2400" "TMe-2271" "TMe-745" 
#> [67] "TMe-1299" "TMe-574"  "TMe-532"  "TMe-2904" "TMe-2769" "TMe-926" 
#> [73] "TMe-340"  "TMe-2057" "TMe-1680" "TMe-332"  "TMe-1880" "TMe-402" 
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-720"  "TMe-726"  "TMe-1481" "TMe-2535" "TMe-1074" "TMe-936" 
#>  [7] "TMe-138"  "TMe-2195" "TMe-1945" "TMe-705"  "TMe-985"  "TMe-1164"
#> [13] "TMe-281"  "TMe-1286" "TMe-983"  "TMe-1554" "TMe-2389" "TMe-505" 
#> [19] "TMe-1847" "TMe-321"  "TMe-852"  "TMe-213"  "TMe-791"  "TMe-1302"
#> [25] "TMe-1518" "TMe-1506" "TMe-1608" "TMe-1744" "TMe-752"  "TMe-662" 
#> [31] "TMe-1704" "TMe-1477" "TMe-1678" "TMe-626"  "TMe-1531" "TMe-1216"
#> [37] "TMe-1832"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_nn_out, use.names = FALSE)) +
  labs(title = "nearest.neighbour")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch selected accessions by cluster based methods
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Globally Optimal Medoid Sampling via Partitioning Around Medoids (PAM)
sel_pam_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "optim.medoid")
sel_pam_out
#> $I
#>  [1] "TMe-1981" "TMe-2810" "TMe-3149" "TMe-132"  "TMe-3424" "TMe-299" 
#>  [7] "TMe-2103" "TMe-2917" "TMe-1086" "TMe-2191" "TMe-489"  "TMe-1716"
#> [13] "TMe-2255" "TMe-3441" "TMe-1091" "TMe-642"  "TMe-736"  "TMe-756" 
#> [19] "TMe-3490" "TMe-28"   "TMe-1823" "TMe-1333" "TMe-1140" "TMe-888" 
#> [25] "TMe-765"  "TMe-1341" "TMe-3202" "TMe-2131" "TMe-2251" "TMe-1532"
#> [31] "TMe-2944" "TMe-910"  "TMe-1621" "TMe-1777" "TMe-1931" "TMe-2066"
#> [37] "TMe-2604" "TMe-3729" "TMe-3425" "TMe-3296" "TMe-3151" "TMe-3220"
#> [43] "TMe-3499" "TMe-3330" "TMe-3170" "TMe-3003" "TMe-3577" "TMe-3685"
#> [49] "TMe-686"  "TMe-2975" "TMe-3229" "TMe-3501" "TMe-2913" "TMe-2936"
#> [55] "TMe-2911" "TMe-3130" "TMe-3419" "TMe-3398" "TMe-3292" "TMe-3389"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2"    "TMe-2611" "TMe-3066" "TMe-85"   "TMe-1505" "TMe-2414"
#>  [7] "TMe-2935" "TMe-681"  "TMe-1242" "TMe-409"  "TMe-3203" "TMe-2242"
#> [13] "TMe-930"  "TMe-2951" "TMe-589"  "TMe-2705" "TMe-2843" "TMe-660" 
#> [19] "TMe-2123" "TMe-3365" "TMe-1833" "TMe-1385" "TMe-1864" "TMe-2048"
#> [25] "TMe-3188" "TMe-1831" "TMe-2077" "TMe-2021" "TMe-2033" "TMe-2568"
#> [31] "TMe-431"  "TMe-74"   "TMe-3093" "TMe-2352" "TMe-3466" "TMe-3605"
#> [37] "TMe-176"  "TMe-1422" "TMe-3210" "TMe-3338" "TMe-3805"
#> 
#> $III
#>  [1] "TMe-3207" "TMe-3028" "TMe-3133" "TMe-1836" "TMe-64"   "TMe-123" 
#>  [7] "TMe-2823" "TMe-3118" "TMe-3731" "TMe-1398" "TMe-206"  "TMe-1138"
#> [13] "TMe-3335" "TMe-2158" "TMe-425"  "TMe-853"  "TMe-2394" "TMe-1790"
#> [19] "TMe-2751" "TMe-3176" "TMe-1725" "TMe-3644" "TMe-3088" "TMe-1897"
#> [25] "TMe-10"   "TMe-1675" "TMe-3643" "TMe-187"  "TMe-3638" "TMe-3324"
#> [31] "TMe-3663" "TMe-3148" "TMe-3598" "TMe-3592" "TMe-3721" "TMe-3166"
#> [37] "TMe-3596"
#> 
#> $IV
#>  [1] "TMe-3073" "TMe-12"   "TMe-204"  "TMe-62"   "TMe-72"   "TMe-525" 
#>  [7] "TMe-1364" "TMe-3440" "TMe-3269" "TMe-207"  "TMe-318"  "TMe-3168"
#> [13] "TMe-1903" "TMe-3107" "TMe-225"  "TMe-513"  "TMe-153"  "TMe-82"  
#> [19] "TMe-380"  "TMe-1150" "TMe-386"  "TMe-387"  "TMe-396"  "TMe-1726"
#> [25] "TMe-2039" "TMe-1765" "TMe-1470" "TMe-3231" "TMe-3435" "TMe-87"  
#> [31] "TMe-550"  "TMe-575"  "TMe-817"  "TMe-1700" "TMe-2410" "TMe-221" 
#> [37] "TMe-1328" "TMe-3541" "TMe-834"  "TMe-190"  "TMe-1336" "TMe-1118"
#> [43] "TMe-1419" "TMe-1139" "TMe-2458" "TMe-802"  "TMe-2567" "TMe-2980"
#> [49] "TMe-2960" "TMe-2210" "TMe-1434" "TMe-3537" "TMe-3357" "TMe-2399"
#> [55] "TMe-2928" "TMe-2947" "TMe-3454" "TMe-78"   "TMe-241"  "TMe-2833"
#> [61] "TMe-2954" "TMe-3585" "TMe-3256" "TMe-3538" "TMe-3198" "TMe-2989"
#> [67] "TMe-3591" "TMe-3068" "TMe-3594" "TMe-3562" "TMe-3273" "TMe-3044"
#> [73] "TMe-3040" "TMe-1278" "TMe-698"  "TMe-699"  "TMe-3253" "TMe-3730"
#> [79] "TMe-3442" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-99"   "TMe-1127" "TMe-1534" "TMe-212"  "TMe-447" 
#>  [7] "TMe-436"  "TMe-479"  "TMe-348"  "TMe-1183" "TMe-582"  "TMe-2009"
#> [13] "TMe-1299" "TMe-3332" "TMe-359"  "TMe-2441" "TMe-826"  "TMe-623" 
#> [19] "TMe-168"  "TMe-378"  "TMe-1011" "TMe-1192" "TMe-417"  "TMe-1234"
#> [25] "TMe-430"  "TMe-1458" "TMe-456"  "TMe-629"  "TMe-1414" "TMe-476" 
#> [31] "TMe-510"  "TMe-2904" "TMe-982"  "TMe-1390" "TMe-667"  "TMe-730" 
#> [37] "TMe-797"  "TMe-669"  "TMe-1901" "TMe-1622" "TMe-1877" "TMe-2213"
#> [43] "TMe-1012" "TMe-700"  "TMe-1760" "TMe-2290" "TMe-2032" "TMe-3411"
#> [49] "TMe-1307" "TMe-1957" "TMe-1440" "TMe-850"  "TMe-1312" "TMe-685" 
#> [55] "TMe-2326" "TMe-1871" "TMe-2192" "TMe-1963" "TMe-2015" "TMe-1487"
#> [61] "TMe-616"  "TMe-2121" "TMe-2435" "TMe-2530" "TMe-1559" "TMe-2790"
#> [67] "TMe-2853" "TMe-217"  "TMe-3185" "TMe-344"  "TMe-536"  "TMe-1227"
#> [73] "TMe-712"  "TMe-771"  "TMe-1733" "TMe-2355" "TMe-2961" "TMe-3408"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-883"  "TMe-2195" "TMe-1832" "TMe-1110" "TMe-852"  "TMe-1512"
#>  [7] "TMe-662"  "TMe-1481" "TMe-2308" "TMe-548"  "TMe-2402" "TMe-1531"
#> [13] "TMe-678"  "TMe-683"  "TMe-1945" "TMe-281"  "TMe-2389" "TMe-1164"
#> [19] "TMe-1477" "TMe-1608" "TMe-2572" "TMe-1017" "TMe-1460" "TMe-625" 
#> [25] "TMe-620"  "TMe-1676" "TMe-1413" "TMe-1416" "TMe-1661" "TMe-1506"
#> [31] "TMe-1992" "TMe-1858" "TMe-1239" "TMe-1319" "TMe-1174" "TMe-1646"
#> [37] "TMe-3095"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_pam_out, use.names = FALSE)) +
  labs(title = "optim.medoid")


# Cluster-Based Sampling via K-means (Naes Method)
sel_naes_out <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "naes")
sel_naes_out
#> $I
#>  [1] "TMe-3418" "TMe-3705" "TMe-2779" "TMe-3461" "TMe-757"  "TMe-2010"
#>  [7] "TMe-306"  "TMe-2131" "TMe-3202" "TMe-248"  "TMe-2191" "TMe-1955"
#> [13] "TMe-410"  "TMe-3394" "TMe-3027" "TMe-3323" "TMe-1621" "TMe-3548"
#> [19] "TMe-3102" "TMe-990"  "TMe-3601" "TMe-3433" "TMe-1981" "TMe-3145"
#> [25] "TMe-3261" "TMe-1190" "TMe-2909" "TMe-2934" "TMe-1096" "TMe-1717"
#> [31] "TMe-3266" "TMe-3416" "TMe-3485" "TMe-1091" "TMe-28"   "TMe-1532"
#> [37] "TMe-3296" "TMe-1156" "TMe-3115" "TMe-501"  "TMe-3625" "TMe-3002"
#> [43] "TMe-3337" "TMe-3208" "TMe-610"  "TMe-2975" "TMe-2943" "TMe-2604"
#> [49] "TMe-3475" "TMe-3425" "TMe-2084" "TMe-3424" "TMe-670"  "TMe-3392"
#> [55] "TMe-2944" "TMe-1930" "TMe-2913" "TMe-1272" "TMe-489"  "TMe-1869"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-160"  "TMe-1385" "TMe-1864" "TMe-2056" "TMe-450"  "TMe-2166"
#>  [7] "TMe-3766" "TMe-589"  "TMe-509"  "TMe-455"  "TMe-1444" "TMe-2970"
#> [13] "TMe-2611" "TMe-1051" "TMe-960"  "TMe-2866" "TMe-3612" "TMe-681" 
#> [19] "TMe-3366" "TMe-674"  "TMe-2258" "TMe-2033" "TMe-1698" "TMe-377" 
#> [25] "TMe-1766" "TMe-2860" "TMe-2824" "TMe-774"  "TMe-3066" "TMe-477" 
#> [31] "TMe-2414" "TMe-3599" "TMe-3052" "TMe-539"  "TMe-2765" "TMe-2304"
#> [37] "TMe-3093" "TMe-1461" "TMe-2951" "TMe-1619" "TMe-2"   
#> 
#> $III
#>  [1] "TMe-3085" "TMe-3336" "TMe-2901" "TMe-13"   "TMe-3032" "TMe-3397"
#>  [7] "TMe-1443" "TMe-3458" "TMe-2481" "TMe-1836" "TMe-1804" "TMe-830" 
#> [13] "TMe-1138" "TMe-3088" "TMe-3148" "TMe-10"   "TMe-187"  "TMe-3324"
#> [19] "TMe-3352" "TMe-946"  "TMe-3796" "TMe-2270" "TMe-3028" "TMe-1262"
#> [25] "TMe-70"   "TMe-853"  "TMe-1910" "TMe-3335" "TMe-1176" "TMe-3721"
#> [31] "TMe-3592" "TMe-3048" "TMe-3544" "TMe-2926" "TMe-2756" "TMe-3804"
#> [37] "TMe-267" 
#> 
#> $IV
#>  [1] "TMe-3752" "TMe-2552" "TMe-3535" "TMe-396"  "TMe-2399" "TMe-513" 
#>  [7] "TMe-2776" "TMe-241"  "TMe-318"  "TMe-314"  "TMe-1765" "TMe-1297"
#> [13] "TMe-1336" "TMe-919"  "TMe-2172" "TMe-3044" "TMe-72"   "TMe-3242"
#> [19] "TMe-2758" "TMe-3660" "TMe-1278" "TMe-1511" "TMe-221"  "TMe-525" 
#> [25] "TMe-1700" "TMe-353"  "TMe-180"  "TMe-3256" "TMe-3454" "TMe-1094"
#> [31] "TMe-3004" "TMe-190"  "TMe-1402" "TMe-634"  "TMe-761"  "TMe-2039"
#> [37] "TMe-3537" "TMe-3167" "TMe-3017" "TMe-3568" "TMe-1419" "TMe-2924"
#> [43] "TMe-3045" "TMe-427"  "TMe-3542" "TMe-885"  "TMe-428"  "TMe-1364"
#> [49] "TMe-204"  "TMe-63"   "TMe-1470" "TMe-2890" "TMe-817"  "TMe-1867"
#> [55] "TMe-3357" "TMe-209"  "TMe-2954" "TMe-1129" "TMe-2775" "TMe-3591"
#> [61] "TMe-1335" "TMe-82"   "TMe-3438" "TMe-2390" "TMe-699"  "TMe-1162"
#> [67] "TMe-283"  "TMe-30"   "TMe-3089" "TMe-1726" "TMe-2210" "TMe-380" 
#> [73] "TMe-1128" "TMe-3780" "TMe-1139" "TMe-1806" "TMe-3297" "TMe-43"  
#> [79] "TMe-3040" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-1572" "TMe-2441" "TMe-323"  "TMe-1358" "TMe-1372" "TMe-804" 
#>  [7] "TMe-1924" "TMe-3363" "TMe-1307" "TMe-284"  "TMe-137"  "TMe-2192"
#> [13] "TMe-424"  "TMe-926"  "TMe-982"  "TMe-224"  "TMe-1629" "TMe-1458"
#> [19] "TMe-167"  "TMe-1188" "TMe-712"  "TMe-823"  "TMe-2769" "TMe-700" 
#> [25] "TMe-2032" "TMe-858"  "TMe-844"  "TMe-1885" "TMe-2213" "TMe-2933"
#> [31] "TMe-447"  "TMe-1715" "TMe-181"  "TMe-1399" "TMe-48"   "TMe-417" 
#> [37] "TMe-418"  "TMe-1778" "TMe-977"  "TMe-294"  "TMe-2820" "TMe-877" 
#> [43] "TMe-343"  "TMe-1418" "TMe-3356" "TMe-870"  "TMe-2296" "TMe-2326"
#> [49] "TMe-439"  "TMe-1322" "TMe-707"  "TMe-1414" "TMe-1339" "TMe-1131"
#> [55] "TMe-1901" "TMe-972"  "TMe-3329" "TMe-861"  "TMe-1227" "TMe-679" 
#> [61] "TMe-543"  "TMe-1199" "TMe-256"  "TMe-1979" "TMe-307"  "TMe-826" 
#> [67] "TMe-1004" "TMe-2853" "TMe-747"  "TMe-2753" "TMe-340"  "TMe-1269"
#> [73] "TMe-142"  "TMe-212"  "TMe-1220" "TMe-1523" "TMe-336"  "TMe-825" 
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1442" "TMe-781"  "TMe-361"  "TMe-3095" "TMe-1506" "TMe-1025"
#>  [7] "TMe-548"  "TMe-1775" "TMe-1796" "TMe-1486" "TMe-1164" "TMe-752" 
#> [13] "TMe-726"  "TMe-2963" "TMe-1124" "TMe-1286" "TMe-852"  "TMe-1646"
#> [19] "TMe-1902" "TMe-2402" "TMe-1769" "TMe-580"  "TMe-2195" "TMe-1477"
#> [25] "TMe-1554" "TMe-269"  "TMe-836"  "TMe-1362" "TMe-838"  "TMe-696" 
#> [31] "TMe-1171" "TMe-659"  "TMe-188"  "TMe-273"  "TMe-1079" "TMe-1558"
#> [37] "TMe-1174"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_naes_out, use.names = FALSE)) +
  labs(title = "naes")


# Cluster-Based Sampling via Hierarchical Clustering with Random Selection

## UPGMA
sel_hclust_random_out1 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random", hclust.method = "average")
sel_hclust_random_out1
#> $I
#>  [1] "TMe-3480" "TMe-41"   "TMe-2810" "TMe-2383" "TMe-2785" "TMe-438" 
#>  [7] "TMe-2069" "TMe-2131" "TMe-410"  "TMe-933"  "TMe-469"  "TMe-715" 
#> [13] "TMe-2027" "TMe-566"  "TMe-896"  "TMe-2779" "TMe-3396" "TMe-2940"
#> [19] "TMe-3521" "TMe-815"  "TMe-867"  "TMe-610"  "TMe-882"  "TMe-839" 
#> [25] "TMe-3461" "TMe-1425" "TMe-2965" "TMe-1621" "TMe-1786" "TMe-2912"
#> [31] "TMe-2010" "TMe-3115" "TMe-2909" "TMe-2967" "TMe-33"   "TMe-3496"
#> [37] "TMe-3577" "TMe-3151" "TMe-3169" "TMe-3076" "TMe-3272" "TMe-3266"
#> [43] "TMe-2916" "TMe-2917" "TMe-3641" "TMe-3229" "TMe-2966" "TMe-2927"
#> [49] "TMe-2984" "TMe-2996" "TMe-3002" "TMe-3112" "TMe-3163" "TMe-3292"
#> [55] "TMe-3389" "TMe-3415" "TMe-3437" "TMe-3471" "TMe-3694" "TMe-3478"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2"    "TMe-706"  "TMe-3371" "TMe-85"   "TMe-3009" "TMe-289" 
#>  [7] "TMe-2978" "TMe-369"  "TMe-2123" "TMe-2085" "TMe-2048" "TMe-2329"
#> [13] "TMe-539"  "TMe-660"  "TMe-902"  "TMe-1005" "TMe-1505" "TMe-1385"
#> [19] "TMe-2412" "TMe-1907" "TMe-1919" "TMe-3805" "TMe-2568" "TMe-2204"
#> [25] "TMe-431"  "TMe-2915" "TMe-3642" "TMe-3020" "TMe-3140" "TMe-3466"
#> [31] "TMe-3605" "TMe-774"  "TMe-1668" "TMe-2952" "TMe-3200" "TMe-3210"
#> [37] "TMe-3222" "TMe-3530" "TMe-3338" "TMe-3406" "TMe-3690"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-1443" "TMe-2977" "TMe-3085" "TMe-3048" "TMe-381" 
#>  [7] "TMe-141"  "TMe-1006" "TMe-3731" "TMe-2169" "TMe-2167" "TMe-261" 
#> [13] "TMe-635"  "TMe-405"  "TMe-2901" "TMe-1897" "TMe-2481" "TMe-946" 
#> [19] "TMe-1200" "TMe-3234" "TMe-3644" "TMe-2270" "TMe-1288" "TMe-2823"
#> [25] "TMe-2897" "TMe-3032" "TMe-3592" "TMe-3397" "TMe-3128" "TMe-3750"
#> [31] "TMe-3299" "TMe-3458" "TMe-3721" "TMe-3608" "TMe-35"   "TMe-2926"
#> [37] "TMe-3143"
#> 
#> $IV
#>  [1] "TMe-3434" "TMe-3231" "TMe-3198" "TMe-7"    "TMe-107"  "TMe-460" 
#>  [7] "TMe-1525" "TMe-2839" "TMe-1027" "TMe-415"  "TMe-2788" "TMe-314" 
#> [13] "TMe-2954" "TMe-956"  "TMe-1267" "TMe-383"  "TMe-1151" "TMe-1010"
#> [19] "TMe-1129" "TMe-394"  "TMe-1456" "TMe-1726" "TMe-427"  "TMe-467" 
#> [25] "TMe-2172" "TMe-513"  "TMe-3435" "TMe-550"  "TMe-3189" "TMe-1652"
#> [31] "TMe-3232" "TMe-1575" "TMe-1406" "TMe-897"  "TMe-919"  "TMe-1095"
#> [37] "TMe-1891" "TMe-2567" "TMe-2980" "TMe-3707" "TMe-186"  "TMe-2956"
#> [43] "TMe-1434" "TMe-3527" "TMe-3242" "TMe-2020" "TMe-2043" "TMe-2240"
#> [49] "TMe-78"   "TMe-2809" "TMe-210"  "TMe-3256" "TMe-3297" "TMe-2946"
#> [55] "TMe-3025" "TMe-3054" "TMe-3053" "TMe-3761" "TMe-3701" "TMe-3273"
#> [61] "TMe-3276" "TMe-3283" "TMe-3428" "TMe-3545" "TMe-3758" "TMe-3773"
#> [67] "TMe-1580" "TMe-3454" "TMe-372"  "TMe-663"  "TMe-698"  "TMe-737" 
#> [73] "TMe-1078" "TMe-1020" "TMe-3097" "TMe-1511" "TMe-3730" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-3659" "TMe-2425" "TMe-1534" "TMe-212"  "TMe-1405"
#>  [7] "TMe-294"  "TMe-323"  "TMe-1643" "TMe-1399" "TMe-969"  "TMe-208" 
#> [13] "TMe-436"  "TMe-826"  "TMe-363"  "TMe-1372" "TMe-976"  "TMe-378" 
#> [19] "TMe-816"  "TMe-584"  "TMe-800"  "TMe-418"  "TMe-667"  "TMe-895" 
#> [25] "TMe-797"  "TMe-932"  "TMe-2307" "TMe-707"  "TMe-592"  "TMe-603" 
#> [31] "TMe-755"  "TMe-819"  "TMe-828"  "TMe-861"  "TMe-1411" "TMe-954" 
#> [37] "TMe-1523" "TMe-1785" "TMe-700"  "TMe-1332" "TMe-1366" "TMe-1440"
#> [43] "TMe-1963" "TMe-1957" "TMe-1629" "TMe-2542" "TMe-1730" "TMe-1541"
#> [49] "TMe-2192" "TMe-1778" "TMe-2041" "TMe-2050" "TMe-2121" "TMe-3034"
#> [55] "TMe-1487" "TMe-2413" "TMe-2590" "TMe-3628" "TMe-2589" "TMe-2790"
#> [61] "TMe-2853" "TMe-3185" "TMe-189"  "TMe-764"  "TMe-344"  "TMe-481" 
#> [67] "TMe-536"  "TMe-712"  "TMe-723"  "TMe-823"  "TMe-1042" "TMe-1196"
#> [73] "TMe-2355" "TMe-2905" "TMe-2961" "TMe-2959" "TMe-3356" "TMe-3736"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-852"  "TMe-1985" "TMe-1076" "TMe-1477" "TMe-963"  "TMe-1403"
#>  [7] "TMe-465"  "TMe-883"  "TMe-620"  "TMe-631"  "TMe-1110" "TMe-1074"
#> [13] "TMe-725"  "TMe-983"  "TMe-1756" "TMe-903"  "TMe-690"  "TMe-1178"
#> [19] "TMe-625"  "TMe-3116" "TMe-1413" "TMe-1416" "TMe-1723" "TMe-1445"
#> [25] "TMe-1506" "TMe-1775" "TMe-1875" "TMe-2035" "TMe-1239" "TMe-2791"
#> [31] "TMe-236"  "TMe-1053" "TMe-1174" "TMe-2551" "TMe-1769" "TMe-2510"
#> [37] "TMe-3095"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out1, use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "average")


## Single-linkage
sel_hclust_random_out2 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random", hclust.method = "single")
sel_hclust_random_out2
#> $I
#>  [1] "TMe-1117" "TMe-41"   "TMe-132"  "TMe-469"  "TMe-500"  "TMe-566" 
#>  [7] "TMe-569"  "TMe-841"  "TMe-815"  "TMe-867"  "TMe-1096" "TMe-1425"
#> [13] "TMe-1466" "TMe-1564" "TMe-1716" "TMe-2027" "TMe-2040" "TMe-2069"
#> [19] "TMe-2372" "TMe-2773" "TMe-2785" "TMe-2967" "TMe-3151" "TMe-3169"
#> [25] "TMe-3262" "TMe-3272" "TMe-3345" "TMe-3353" "TMe-3452" "TMe-3633"
#> [31] "TMe-606"  "TMe-715"  "TMe-717"  "TMe-1218" "TMe-2462" "TMe-2912"
#> [37] "TMe-2934" "TMe-2937" "TMe-2984" "TMe-2985" "TMe-2993" "TMe-3112"
#> [43] "TMe-3163" "TMe-3229" "TMe-3282" "TMe-3292" "TMe-3389" "TMe-3415"
#> [49] "TMe-3433" "TMe-3437" "TMe-3471" "TMe-3475" "TMe-3488" "TMe-3514"
#> [55] "TMe-3694" "TMe-2940" "TMe-2943" "TMe-3462" "TMe-3463" "TMe-3478"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3354" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-404"  "TMe-768" 
#>  [7] "TMe-902"  "TMe-1242" "TMe-1385" "TMe-1444" "TMe-1831" "TMe-1907"
#> [13] "TMe-1919" "TMe-2033" "TMe-2054" "TMe-2056" "TMe-2204" "TMe-2851"
#> [19] "TMe-2866" "TMe-2891" "TMe-3009" "TMe-3140" "TMe-3766" "TMe-3466"
#> [25] "TMe-3605" "TMe-47"   "TMe-176"  "TMe-431"  "TMe-478"  "TMe-648" 
#> [31] "TMe-1310" "TMe-2352" "TMe-2568" "TMe-2952" "TMe-3200" "TMe-3210"
#> [37] "TMe-3222" "TMe-3338" "TMe-3690" "TMe-3771" "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-3007" "TMe-104"  "TMe-261"  "TMe-267" 
#>  [7] "TMe-381"  "TMe-425"  "TMe-785"  "TMe-1200" "TMe-1593" "TMe-1897"
#> [13] "TMe-2270" "TMe-2374" "TMe-2823" "TMe-2897" "TMe-3299" "TMe-3804"
#> [19] "TMe-35"   "TMe-116"  "TMe-174"  "TMe-187"  "TMe-234"  "TMe-635" 
#> [25] "TMe-2756" "TMe-3008" "TMe-3143" "TMe-3274" "TMe-3383" "TMe-3397"
#> [31] "TMe-3445" "TMe-3551" "TMe-3565" "TMe-3575" "TMe-3596" "TMe-3620"
#> [37] "TMe-3721"
#> 
#> $IV
#>  [1] "TMe-3000" "TMe-154"  "TMe-314"  "TMe-427"  "TMe-650"  "TMe-766" 
#>  [7] "TMe-812"  "TMe-897"  "TMe-919"  "TMe-956"  "TMe-962"  "TMe-967" 
#> [13] "TMe-1342" "TMe-1350" "TMe-1368" "TMe-1369" "TMe-1376" "TMe-1406"
#> [19] "TMe-1434" "TMe-1525" "TMe-1806" "TMe-2020" "TMe-2025" "TMe-2043"
#> [25] "TMe-241"  "TMe-2839" "TMe-2956" "TMe-2971" "TMe-3025" "TMe-3054"
#> [31] "TMe-3055" "TMe-3089" "TMe-3255" "TMe-3256" "TMe-3273" "TMe-3275"
#> [37] "TMe-3591" "TMe-3545" "TMe-5"    "TMe-25"   "TMe-180"  "TMe-205" 
#> [43] "TMe-209"  "TMe-225"  "TMe-226"  "TMe-259"  "TMe-280"  "TMe-372" 
#> [49] "TMe-608"  "TMe-619"  "TMe-663"  "TMe-698"  "TMe-733"  "TMe-735" 
#> [55] "TMe-737"  "TMe-1010" "TMe-1020" "TMe-1129" "TMe-3109" "TMe-1181"
#> [61] "TMe-1511" "TMe-1988" "TMe-2567" "TMe-2922" "TMe-3167" "TMe-3189"
#> [67] "TMe-3225" "TMe-3228" "TMe-3297" "TMe-3386" "TMe-3422" "TMe-3527"
#> [73] "TMe-3701" "TMe-3707" "TMe-3730" "TMe-3053" "TMe-3316" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-162"  "TMe-347"  "TMe-376"  "TMe-378"  "TMe-399" 
#>  [7] "TMe-418"  "TMe-419"  "TMe-423"  "TMe-782"  "TMe-819"  "TMe-861" 
#> [13] "TMe-972"  "TMe-1099" "TMe-1204" "TMe-1357" "TMe-1366" "TMe-1381"
#> [19] "TMe-1458" "TMe-1482" "TMe-1488" "TMe-1730" "TMe-1873" "TMe-1924"
#> [25] "TMe-1963" "TMe-2041" "TMe-2050" "TMe-2060" "TMe-2121" "TMe-2192"
#> [31] "TMe-2290" "TMe-2413" "TMe-2531" "TMe-2589" "TMe-2612" "TMe-2688"
#> [37] "TMe-2750" "TMe-2790" "TMe-2853" "TMe-2863" "TMe-2973" "TMe-3185"
#> [43] "TMe-189"  "TMe-194"  "TMe-208"  "TMe-277"  "TMe-344"  "TMe-341" 
#> [49] "TMe-481"  "TMe-536"  "TMe-616"  "TMe-700"  "TMe-712"  "TMe-723" 
#> [55] "TMe-832"  "TMe-866"  "TMe-877"  "TMe-1042" "TMe-1196" "TMe-1199"
#> [61] "TMe-1300" "TMe-1305" "TMe-1307" "TMe-1311" "TMe-1391" "TMe-1453"
#> [67] "TMe-2003" "TMe-2139" "TMe-2355" "TMe-2905" "TMe-2953" "TMe-2959"
#> [73] "TMe-3034" "TMe-3311" "TMe-3356" "TMe-3363" "TMe-3408" "TMe-3736"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-130"  "TMe-1985" "TMe-269"  "TMe-983"  "TMe-361"  "TMe-442" 
#>  [7] "TMe-465"  "TMe-514"  "TMe-1383" "TMe-1416" "TMe-1428" "TMe-1518"
#> [13] "TMe-1775" "TMe-1875" "TMe-2035" "TMe-2064" "TMe-2067" "TMe-2196"
#> [19] "TMe-2534" "TMe-690"  "TMe-693"  "TMe-1053" "TMe-1124" "TMe-1174"
#> [25] "TMe-1182" "TMe-1232" "TMe-1239" "TMe-2551" "TMe-1723" "TMe-1769"
#> [31] "TMe-2510" "TMe-2963" "TMe-3095" "TMe-3387" "TMe-2983" "TMe-3116"
#> [37] "TMe-3549"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out2, use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "single")


## Complete-linkage
sel_hclust_random_out3 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random", hclust.method = "complete")
sel_hclust_random_out3
#> $I
#>  [1] "TMe-642"  "TMe-1621" "TMe-469"  "TMe-1600" "TMe-1830" "TMe-438" 
#>  [7] "TMe-1981" "TMe-1973" "TMe-410"  "TMe-1589" "TMe-1156" "TMe-500" 
#> [13] "TMe-2226" "TMe-1360" "TMe-1958" "TMe-3496" "TMe-736"  "TMe-3249"
#> [19] "TMe-1716" "TMe-778"  "TMe-3345" "TMe-3264" "TMe-1533" "TMe-940" 
#> [25] "TMe-1960" "TMe-579"  "TMe-1272" "TMe-1117" "TMe-2372" "TMe-2191"
#> [31] "TMe-1914" "TMe-2912" "TMe-1532" "TMe-3319" "TMe-2773" "TMe-3518"
#> [37] "TMe-2984" "TMe-3065" "TMe-3030" "TMe-2914" "TMe-3416" "TMe-3110"
#> [43] "TMe-3111" "TMe-3479" "TMe-3433" "TMe-3220" "TMe-2985" "TMe-3330"
#> [49] "TMe-2909" "TMe-3625" "TMe-3641" "TMe-686"  "TMe-2910" "TMe-3471"
#> [55] "TMe-3292" "TMe-3102" "TMe-3389" "TMe-3461" "TMe-3724" "TMe-3601"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3150" "TMe-3401" "TMe-3066" "TMe-3466" "TMe-251"  "TMe-1732"
#>  [7] "TMe-3239" "TMe-478"  "TMe-1323" "TMe-2428" "TMe-664"  "TMe-890" 
#> [13] "TMe-477"  "TMe-589"  "TMe-2951" "TMe-1504" "TMe-2217" "TMe-2268"
#> [19] "TMe-3365" "TMe-1754" "TMe-1005" "TMe-2123" "TMe-3690" "TMe-1385"
#> [25] "TMe-3146" "TMe-2021" "TMe-1271" "TMe-2033" "TMe-2568" "TMe-2995"
#> [31] "TMe-1459" "TMe-2915" "TMe-3642" "TMe-3093" "TMe-3020" "TMe-2352"
#> [37] "TMe-2860" "TMe-1616" "TMe-2952" "TMe-3366" "TMe-3530"
#> 
#> $III
#>  [1] "TMe-155"  "TMe-200"  "TMe-14"   "TMe-64"   "TMe-206"  "TMe-123" 
#>  [7] "TMe-2751" "TMe-1956" "TMe-2154" "TMe-1725" "TMe-635"  "TMe-420" 
#> [13] "TMe-1868" "TMe-1738" "TMe-2405" "TMe-1787" "TMe-1176" "TMe-2756"
#> [19] "TMe-3302" "TMe-1863" "TMe-2811" "TMe-3007" "TMe-3608" "TMe-3324"
#> [25] "TMe-2968" "TMe-3663" "TMe-2926" "TMe-3010" "TMe-3383" "TMe-3299"
#> [31] "TMe-3458" "TMe-3804" "TMe-23"   "TMe-234"  "TMe-3575" "TMe-3445"
#> [37] "TMe-3565"
#> 
#> $IV
#>  [1] "TMe-3535" "TMe-12"   "TMe-204"  "TMe-3773" "TMe-103"  "TMe-1094"
#>  [7] "TMe-298"  "TMe-154"  "TMe-54"   "TMe-207"  "TMe-235"  "TMe-270" 
#> [13] "TMe-65"   "TMe-225"  "TMe-317"  "TMe-3228" "TMe-415"  "TMe-368" 
#> [19] "TMe-380"  "TMe-1166" "TMe-2890" "TMe-444"  "TMe-3109" "TMe-3639"
#> [25] "TMe-734"  "TMe-609"  "TMe-1765" "TMe-556"  "TMe-2039" "TMe-526" 
#> [31] "TMe-540"  "TMe-184"  "TMe-3542" "TMe-1553" "TMe-1479" "TMe-2172"
#> [37] "TMe-675"  "TMe-2946" "TMe-3573" "TMe-3541" "TMe-2788" "TMe-1278"
#> [43] "TMe-1988" "TMe-699"  "TMe-2399" "TMe-2960" "TMe-3438" "TMe-1434"
#> [49] "TMe-1987" "TMe-1776" "TMe-3218" "TMe-2025" "TMe-1916" "TMe-3545"
#> [55] "TMe-2043" "TMe-2059" "TMe-3537" "TMe-241"  "TMe-737"  "TMe-24"  
#> [61] "TMe-3068" "TMe-3484" "TMe-3730" "TMe-3055" "TMe-3089" "TMe-3255"
#> [67] "TMe-3275" "TMe-3283" "TMe-3615" "TMe-27"   "TMe-186"  "TMe-15"  
#> [73] "TMe-3454" "TMe-684"  "TMe-698"  "TMe-1020" "TMe-3097" "TMe-3707"
#> [79] "TMe-3417" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-98"   "TMe-158"  "TMe-167"  "TMe-355"  "TMe-326" 
#>  [7] "TMe-2750" "TMe-479"  "TMe-3659" "TMe-473"  "TMe-712"  "TMe-2339"
#> [13] "TMe-2643" "TMe-436"  "TMe-330"  "TMe-1877" "TMe-866"  "TMe-2959"
#> [19] "TMe-795"  "TMe-391"  "TMe-901"  "TMe-399"  "TMe-823"  "TMe-2863"
#> [25] "TMe-759"  "TMe-822"  "TMe-667"  "TMe-1859" "TMe-1834" "TMe-870" 
#> [31] "TMe-1522" "TMe-1359" "TMe-543"  "TMe-1234" "TMe-929"  "TMe-755" 
#> [37] "TMe-826"  "TMe-687"  "TMe-1643" "TMe-2355" "TMe-977"  "TMe-287" 
#> [43] "TMe-1357" "TMe-3329" "TMe-1381" "TMe-162"  "TMe-1979" "TMe-1482"
#> [49] "TMe-1957" "TMe-926"  "TMe-2688" "TMe-1730" "TMe-1871" "TMe-2307"
#> [55] "TMe-1311" "TMe-2050" "TMe-2121" "TMe-2151" "TMe-3034" "TMe-2003"
#> [61] "TMe-2435" "TMe-2530" "TMe-2589" "TMe-2853" "TMe-3185" "TMe-1722"
#> [67] "TMe-189"  "TMe-217"  "TMe-344"  "TMe-800"  "TMe-341"  "TMe-397" 
#> [73] "TMe-1199" "TMe-536"  "TMe-723"  "TMe-1778" "TMe-1785" "TMe-3311"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-130"  "TMe-1985" "TMe-1017" "TMe-1217" "TMe-725"  "TMe-442" 
#>  [7] "TMe-465"  "TMe-1531" "TMe-1809" "TMe-661"  "TMe-580"  "TMe-598" 
#> [13] "TMe-720"  "TMe-751"  "TMe-856"  "TMe-846"  "TMe-1238" "TMe-1404"
#> [19] "TMe-1264" "TMe-1261" "TMe-1900" "TMe-1362" "TMe-3095" "TMe-1483"
#> [25] "TMe-1416" "TMe-1661" "TMe-693"  "TMe-2067" "TMe-1506" "TMe-2510"
#> [31] "TMe-1614" "TMe-2983" "TMe-1608" "TMe-236"  "TMe-1053" "TMe-1174"
#> [37] "TMe-2551"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out3, use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "complete")


## Ward's D
sel_hclust_random_out4 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random", hclust.method = "ward.D")
sel_hclust_random_out4
#> $I
#>  [1] "TMe-841"  "TMe-2773" "TMe-3149" "TMe-132"  "TMe-717"  "TMe-2453"
#>  [7] "TMe-1564" "TMe-489"  "TMe-3202" "TMe-1716" "TMe-3496" "TMe-3441"
#> [13] "TMe-569"  "TMe-3480" "TMe-736"  "TMe-3396" "TMe-757"  "TMe-778" 
#> [19] "TMe-3164" "TMe-3266" "TMe-1935" "TMe-610"  "TMe-1096" "TMe-1272"
#> [25] "TMe-1091" "TMe-579"  "TMe-2237" "TMe-1425" "TMe-2251" "TMe-1117"
#> [31] "TMe-1621" "TMe-1869" "TMe-2912" "TMe-3359" "TMe-2152" "TMe-3319"
#> [37] "TMe-3074" "TMe-3211" "TMe-3729" "TMe-3169" "TMe-3111" "TMe-3481"
#> [43] "TMe-3220" "TMe-3235" "TMe-2985" "TMe-2914" "TMe-3488" "TMe-3229"
#> [49] "TMe-3685" "TMe-686"  "TMe-3292" "TMe-3501" "TMe-2927" "TMe-3493"
#> [55] "TMe-3394" "TMe-3694" "TMe-3076" "TMe-3475" "TMe-3548" "TMe-3281"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-40"   "TMe-3284" "TMe-3605" "TMe-85"   "TMe-3009" "TMe-289" 
#>  [7] "TMe-2935" "TMe-681"  "TMe-478"  "TMe-2705" "TMe-1271" "TMe-2085"
#> [13] "TMe-3800" "TMe-1444" "TMe-528"  "TMe-3258" "TMe-589"  "TMe-1504"
#> [19] "TMe-636"  "TMe-2268" "TMe-1864" "TMe-1005" "TMe-1241" "TMe-1172"
#> [25] "TMe-2997" "TMe-171"  "TMe-2242" "TMe-2166" "TMe-2412" "TMe-2033"
#> [31] "TMe-2056" "TMe-1459" "TMe-2851" "TMe-1310" "TMe-3642" "TMe-2978"
#> [37] "TMe-3690" "TMe-3406" "TMe-3141" "TMe-3338" "TMe-3188"
#> 
#> $III
#>  [1] "TMe-2161" "TMe-6"    "TMe-3644" "TMe-304"  "TMe-64"   "TMe-3772"
#>  [7] "TMe-1006" "TMe-405"  "TMe-206"  "TMe-1797" "TMe-261"  "TMe-635" 
#> [13] "TMe-420"  "TMe-1897" "TMe-773"  "TMe-2347" "TMe-1868" "TMe-1593"
#> [19] "TMe-1939" "TMe-2733" "TMe-3277" "TMe-3287" "TMe-1421" "TMe-2356"
#> [25] "TMe-3094" "TMe-2748" "TMe-35"   "TMe-3274" "TMe-3043" "TMe-3230"
#> [31] "TMe-3128" "TMe-3638" "TMe-3598" "TMe-3166" "TMe-3143" "TMe-3445"
#> [37] "TMe-3721"
#> 
#> $IV
#>  [1] "TMe-11"   "TMe-3278" "TMe-3198" "TMe-1765" "TMe-81"   "TMe-806" 
#>  [7] "TMe-1525" "TMe-1577" "TMe-182"  "TMe-318"  "TMe-415"  "TMe-215" 
#> [13] "TMe-735"  "TMe-3760" "TMe-1406" "TMe-386"  "TMe-1181" "TMe-3109"
#> [19] "TMe-2378" "TMe-403"  "TMe-427"  "TMe-3214" "TMe-2172" "TMe-3390"
#> [25] "TMe-3435" "TMe-87"   "TMe-550"  "TMe-3189" "TMe-641"  "TMe-2788"
#> [31] "TMe-1144" "TMe-760"  "TMe-372"  "TMe-1368" "TMe-210"  "TMe-3422"
#> [37] "TMe-1139" "TMe-1820" "TMe-2039" "TMe-320"  "TMe-698"  "TMe-1419"
#> [43] "TMe-3752" "TMe-1297" "TMe-1342" "TMe-2980" "TMe-3707" "TMe-733" 
#> [49] "TMe-2956" "TMe-1374" "TMe-1179" "TMe-1806" "TMe-3017" "TMe-3537"
#> [55] "TMe-38"   "TMe-1903" "TMe-2979" "TMe-2043" "TMe-1397" "TMe-2958"
#> [61] "TMe-241"  "TMe-2988" "TMe-3781" "TMe-3025" "TMe-3103" "TMe-3072"
#> [67] "TMe-3594" "TMe-2946" "TMe-3107" "TMe-3562" "TMe-3290" "TMe-3701"
#> [73] "TMe-3259" "TMe-3773" "TMe-3312" "TMe-15"   "TMe-2758" "TMe-714" 
#> [79] "TMe-3417" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-1196" "TMe-3659" "TMe-158"  "TMe-167"  "TMe-212"  "TMe-868" 
#>  [7] "TMe-294"  "TMe-479"  "TMe-340"  "TMe-334"  "TMe-951"  "TMe-208" 
#> [13] "TMe-197"  "TMe-2361" "TMe-2518" "TMe-1401" "TMe-373"  "TMe-804" 
#> [19] "TMe-2959" "TMe-665"  "TMe-391"  "TMe-390"  "TMe-1880" "TMe-800" 
#> [25] "TMe-406"  "TMe-1042" "TMe-217"  "TMe-2121" "TMe-2030" "TMe-424" 
#> [31] "TMe-832"  "TMe-889"  "TMe-1715" "TMe-476"  "TMe-1009" "TMe-870" 
#> [37] "TMe-543"  "TMe-982"  "TMe-565"  "TMe-592"  "TMe-667"  "TMe-730" 
#> [43] "TMe-750"  "TMe-2853" "TMe-1390" "TMe-1643" "TMe-1877" "TMe-977" 
#> [49] "TMe-2790" "TMe-813"  "TMe-1204" "TMe-2290" "TMe-168"  "TMe-1357"
#> [55] "TMe-1366" "TMe-1375" "TMe-1440" "TMe-1488" "TMe-1694" "TMe-742" 
#> [61] "TMe-2835" "TMe-1285" "TMe-1871" "TMe-1391" "TMe-1934" "TMe-1311"
#> [67] "TMe-2050" "TMe-850"  "TMe-798"  "TMe-2296" "TMe-194"  "TMe-3185"
#> [73] "TMe-245"  "TMe-344"  "TMe-1227" "TMe-747"  "TMe-1778" "TMe-2905"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-130"  "TMe-2195" "TMe-310"  "TMe-809"  "TMe-361"  "TMe-1809"
#>  [7] "TMe-2308" "TMe-470"  "TMe-531"  "TMe-1110" "TMe-1883" "TMe-185" 
#> [13] "TMe-751"  "TMe-2535" "TMe-985"  "TMe-838"  "TMe-1756" "TMe-903" 
#> [19] "TMe-968"  "TMe-1460" "TMe-690"  "TMe-2026" "TMe-1592" "TMe-1124"
#> [25] "TMe-1661" "TMe-1614" "TMe-1506" "TMe-1566" "TMe-1775" "TMe-1835"
#> [31] "TMe-1875" "TMe-2983" "TMe-2196" "TMe-1319" "TMe-1053" "TMe-1608"
#> [37] "TMe-3095"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out4, use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "ward.D")


## WPGMA
sel_hclust_random_out5 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random", hclust.method = "mcquitty")
sel_hclust_random_out5
#> $I
#>  [1] "TMe-1981" "TMe-2810" "TMe-117"  "TMe-132"  "TMe-765"  "TMe-300" 
#>  [7] "TMe-1890" "TMe-410"  "TMe-1906" "TMe-3202" "TMe-1096" "TMe-500" 
#> [13] "TMe-3030" "TMe-566"  "TMe-2191" "TMe-642"  "TMe-3427" "TMe-3314"
#> [19] "TMe-757"  "TMe-778"  "TMe-815"  "TMe-867"  "TMe-940"  "TMe-670" 
#> [25] "TMe-3441" "TMe-1423" "TMe-1964" "TMe-3719" "TMe-1786" "TMe-3729"
#> [31] "TMe-2912" "TMe-3319" "TMe-3457" "TMe-3521" "TMe-2967" "TMe-3515"
#> [37] "TMe-3151" "TMe-3601" "TMe-3577" "TMe-3272" "TMe-3266" "TMe-3465"
#> [43] "TMe-3685" "TMe-717"  "TMe-1532" "TMe-2462" "TMe-3485" "TMe-2916"
#> [49] "TMe-2927" "TMe-33"   "TMe-2984" "TMe-3026" "TMe-3163" "TMe-3292"
#> [55] "TMe-3389" "TMe-3415" "TMe-3437" "TMe-3471" "TMe-3694" "TMe-2910"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3540" "TMe-2211" "TMe-3617" "TMe-3009" "TMe-289"  "TMe-2851"
#>  [7] "TMe-369"  "TMe-1242" "TMe-1754" "TMe-2242" "TMe-477"  "TMe-589" 
#> [13] "TMe-539"  "TMe-3209" "TMe-2268" "TMe-2797" "TMe-902"  "TMe-1005"
#> [19] "TMe-2414" "TMe-1187" "TMe-1385" "TMe-2428" "TMe-1827" "TMe-3447"
#> [25] "TMe-1860" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2056" "TMe-1459"
#> [31] "TMe-3531" "TMe-3766" "TMe-3466" "TMe-3605" "TMe-3366" "TMe-3200"
#> [37] "TMe-3210" "TMe-3530" "TMe-3338" "TMe-3690" "TMe-3800"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-3133" "TMe-3700" "TMe-3544" "TMe-161" 
#>  [7] "TMe-70"   "TMe-1819" "TMe-187"  "TMe-155"  "TMe-3287" "TMe-2405"
#> [13] "TMe-261"  "TMe-267"  "TMe-405"  "TMe-420"  "TMe-1897" "TMe-1738"
#> [19] "TMe-785"  "TMe-2756" "TMe-1790" "TMe-2592" "TMe-3804" "TMe-1939"
#> [25] "TMe-2362" "TMe-2897" "TMe-3382" "TMe-3032" "TMe-3230" "TMe-3565"
#> [31] "TMe-3376" "TMe-3176" "TMe-3299" "TMe-3638" "TMe-3608" "TMe-2926"
#> [37] "TMe-3143"
#> 
#> $IV
#>  [1] "TMe-3206" "TMe-18"   "TMe-3590" "TMe-3243" "TMe-1765" "TMe-54"  
#>  [7] "TMe-534"  "TMe-1246" "TMe-154"  "TMe-1027" "TMe-270"  "TMe-314" 
#> [13] "TMe-3228" "TMe-956"  "TMe-1456" "TMe-1814" "TMe-2399" "TMe-1010"
#> [19] "TMe-766"  "TMe-1148" "TMe-3511" "TMe-1726" "TMe-427"  "TMe-1975"
#> [25] "TMe-1123" "TMe-540"  "TMe-550"  "TMe-1479" "TMe-699"  "TMe-812" 
#> [31] "TMe-897"  "TMe-3422" "TMe-619"  "TMe-1923" "TMe-1162" "TMe-1806"
#> [37] "TMe-2318" "TMe-608"  "TMe-1348" "TMe-3541" "TMe-1351" "TMe-2960"
#> [43] "TMe-885"  "TMe-1434" "TMe-1480" "TMe-1489" "TMe-226"  "TMe-38"  
#> [49] "TMe-2337" "TMe-3734" "TMe-3545" "TMe-2043" "TMe-2958" "TMe-5"   
#> [55] "TMe-2809" "TMe-2833" "TMe-3538" "TMe-3108" "TMe-3025" "TMe-3054"
#> [61] "TMe-3055" "TMe-3761" "TMe-3535" "TMe-3255" "TMe-3273" "TMe-3253"
#> [67] "TMe-3591" "TMe-3581" "TMe-2758" "TMe-3602" "TMe-698"  "TMe-737" 
#> [73] "TMe-761"  "TMe-1020" "TMe-1597" "TMe-3189" "TMe-3730" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-98"   "TMe-2425" "TMe-1534" "TMe-647"  "TMe-651" 
#>  [7] "TMe-2750" "TMe-975"  "TMe-359"  "TMe-1283" "TMe-334"  "TMe-336" 
#> [13] "TMe-193"  "TMe-362"  "TMe-1382" "TMe-1388" "TMe-378"  "TMe-629" 
#> [19] "TMe-391"  "TMe-390"  "TMe-816"  "TMe-418"  "TMe-1458" "TMe-1257"
#> [25] "TMe-797"  "TMe-1440" "TMe-2907" "TMe-1199" "TMe-574"  "TMe-588" 
#> [31] "TMe-803"  "TMe-603"  "TMe-755"  "TMe-2151" "TMe-612"  "TMe-861" 
#> [37] "TMe-920"  "TMe-997"  "TMe-977"  "TMe-972"  "TMe-700"  "TMe-1192"
#> [43] "TMe-1523" "TMe-3411" "TMe-1307" "TMe-1427" "TMe-1963" "TMe-1488"
#> [49] "TMe-2688" "TMe-2542" "TMe-2753" "TMe-1871" "TMe-1934" "TMe-1269"
#> [55] "TMe-2050" "TMe-2121" "TMe-2290" "TMe-2296" "TMe-2413" "TMe-2820"
#> [61] "TMe-2531" "TMe-2589" "TMe-2853" "TMe-2973" "TMe-3185" "TMe-771" 
#> [67] "TMe-284"  "TMe-832"  "TMe-2953" "TMe-712"  "TMe-723"  "TMe-1042"
#> [73] "TMe-1785" "TMe-2905" "TMe-2961" "TMe-2959" "TMe-3356" "TMe-3736"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-130"  "TMe-2195" "TMe-1264" "TMe-470"  "TMe-1503" "TMe-442" 
#>  [7] "TMe-1232" "TMe-2308" "TMe-548"  "TMe-1110" "TMe-2402" "TMe-968" 
#> [13] "TMe-725"  "TMe-1025" "TMe-222"  "TMe-845"  "TMe-1079" "TMe-1472"
#> [19] "TMe-1261" "TMe-1353" "TMe-1383" "TMe-1483" "TMe-1124" "TMe-1445"
#> [25] "TMe-1518" "TMe-1566" "TMe-1900" "TMe-2551" "TMe-2035" "TMe-1744"
#> [31] "TMe-2196" "TMe-2791" "TMe-1319" "TMe-1053" "TMe-1174" "TMe-1769"
#> [37] "TMe-3549"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out5, use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "mcquitty")


## WPGMC
sel_hclust_random_out6 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random", hclust.method = "median")
sel_hclust_random_out6
#> $I
#>  [1] "TMe-2462" "TMe-41"   "TMe-44"   "TMe-117"  "TMe-410"  "TMe-469" 
#>  [7] "TMe-500"  "TMe-566"  "TMe-569"  "TMe-3264" "TMe-815"  "TMe-1096"
#> [13] "TMe-765"  "TMe-1425" "TMe-1451" "TMe-1466" "TMe-1564" "TMe-1922"
#> [19] "TMe-1935" "TMe-1958" "TMe-1973" "TMe-2010" "TMe-2027" "TMe-2372"
#> [25] "TMe-2773" "TMe-2967" "TMe-3151" "TMe-3272" "TMe-3345" "TMe-3353"
#> [31] "TMe-3452" "TMe-3481" "TMe-3633" "TMe-606"  "TMe-717"  "TMe-1218"
#> [37] "TMe-1589" "TMe-2912" "TMe-2927" "TMe-2934" "TMe-2937" "TMe-2985"
#> [43] "TMe-2996" "TMe-3026" "TMe-3112" "TMe-3163" "TMe-3292" "TMe-3389"
#> [49] "TMe-3415" "TMe-3433" "TMe-3437" "TMe-3457" "TMe-3471" "TMe-3475"
#> [55] "TMe-3480" "TMe-3514" "TMe-3705" "TMe-2910" "TMe-2943" "TMe-3462"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2903" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-902"  "TMe-509" 
#>  [7] "TMe-768"  "TMe-1444" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2568"
#> [13] "TMe-2128" "TMe-2204" "TMe-2891" "TMe-2970" "TMe-3009" "TMe-3140"
#> [19] "TMe-3766" "TMe-3401" "TMe-3466" "TMe-3605" "TMe-47"   "TMe-1619"
#> [25] "TMe-2352" "TMe-2860" "TMe-2952" "TMe-2957" "TMe-3200" "TMe-3210"
#> [31] "TMe-3222" "TMe-3237" "TMe-3338" "TMe-3366" "TMe-3406" "TMe-3530"
#> [37] "TMe-3547" "TMe-3690" "TMe-3800" "TMe-3286" "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-2901" "TMe-104"  "TMe-141"  "TMe-261" 
#>  [7] "TMe-267"  "TMe-381"  "TMe-405"  "TMe-420"  "TMe-425"  "TMe-785" 
#> [13] "TMe-946"  "TMe-1200" "TMe-1868" "TMe-1897" "TMe-1993" "TMe-2802"
#> [19] "TMe-2823" "TMe-2897" "TMe-2968" "TMe-3234" "TMe-3299" "TMe-3804"
#> [25] "TMe-35"   "TMe-174"  "TMe-234"  "TMe-635"  "TMe-2756" "TMe-3596"
#> [31] "TMe-3143" "TMe-3274" "TMe-3383" "TMe-3397" "TMe-3445" "TMe-3565"
#> [37] "TMe-3575"
#> 
#> $IV
#>  [1] "TMe-108"  "TMe-154"  "TMe-266"  "TMe-2960" "TMe-427"  "TMe-540" 
#>  [7] "TMe-550"  "TMe-656"  "TMe-734"  "TMe-812"  "TMe-897"  "TMe-919" 
#> [13] "TMe-1162" "TMe-1350" "TMe-1376" "TMe-1406" "TMe-1434" "TMe-1525"
#> [19] "TMe-1776" "TMe-1806" "TMe-1820" "TMe-2025" "TMe-2043" "TMe-2755"
#> [25] "TMe-241"  "TMe-2956" "TMe-2971" "TMe-2979" "TMe-3025" "TMe-3054"
#> [31] "TMe-3055" "TMe-3089" "TMe-3144" "TMe-3196" "TMe-3255" "TMe-3256"
#> [37] "TMe-3273" "TMe-3428" "TMe-3545" "TMe-3576" "TMe-3752" "TMe-5"   
#> [43] "TMe-27"   "TMe-205"  "TMe-209"  "TMe-225"  "TMe-226"  "TMe-259" 
#> [49] "TMe-372"  "TMe-398"  "TMe-421"  "TMe-608"  "TMe-619"  "TMe-663" 
#> [55] "TMe-698"  "TMe-737"  "TMe-761"  "TMe-1010" "TMe-1020" "TMe-1078"
#> [61] "TMe-1129" "TMe-1511" "TMe-1597" "TMe-1988" "TMe-3097" "TMe-3189"
#> [67] "TMe-3297" "TMe-3386" "TMe-3422" "TMe-3527" "TMe-3558" "TMe-3619"
#> [73] "TMe-3701" "TMe-3707" "TMe-3727" "TMe-3730" "TMe-3053" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-2612" "TMe-2750" "TMe-969"  "TMe-362"  "TMe-373" 
#>  [7] "TMe-376"  "TMe-378"  "TMe-543"  "TMe-399"  "TMe-418"  "TMe-419" 
#> [13] "TMe-423"  "TMe-547"  "TMe-603"  "TMe-755"  "TMe-730"  "TMe-828" 
#> [19] "TMe-861"  "TMe-1099" "TMe-1204" "TMe-1366" "TMe-1375" "TMe-1381"
#> [25] "TMe-1427" "TMe-1458" "TMe-1482" "TMe-1963" "TMe-2041" "TMe-2050"
#> [31] "TMe-2060" "TMe-2121" "TMe-2531" "TMe-2589" "TMe-2688" "TMe-2790"
#> [37] "TMe-2853" "TMe-2855" "TMe-2863" "TMe-2973" "TMe-3185" "TMe-162" 
#> [43] "TMe-189"  "TMe-193"  "TMe-208"  "TMe-224"  "TMe-277"  "TMe-481" 
#> [49] "TMe-536"  "TMe-616"  "TMe-627"  "TMe-688"  "TMe-700"  "TMe-712" 
#> [55] "TMe-723"  "TMe-748"  "TMe-816"  "TMe-866"  "TMe-975"  "TMe-1042"
#> [61] "TMe-1196" "TMe-1283" "TMe-1332" "TMe-1453" "TMe-1785" "TMe-2003"
#> [67] "TMe-2139" "TMe-2355" "TMe-2905" "TMe-2953" "TMe-2959" "TMe-3034"
#> [73] "TMe-3311" "TMe-3356" "TMe-3408" "TMe-3628" "TMe-3736" "TMe-3659"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-854"  "TMe-269"  "TMe-315"  "TMe-442"  "TMe-531"  "TMe-580" 
#>  [7] "TMe-725"  "TMe-751"  "TMe-1383" "TMe-1403" "TMe-1416" "TMe-1428"
#> [13] "TMe-1445" "TMe-1518" "TMe-1566" "TMe-1775" "TMe-1875" "TMe-1995"
#> [19] "TMe-2026" "TMe-2035" "TMe-2067" "TMe-2534" "TMe-2791" "TMe-310" 
#> [25] "TMe-693"  "TMe-963"  "TMe-1053" "TMe-1174" "TMe-1232" "TMe-1558"
#> [31] "TMe-1560" "TMe-1769" "TMe-2510" "TMe-2543" "TMe-3549" "TMe-2983"
#> [37] "TMe-3116"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out6, use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "median")


## UPGMC
sel_hclust_random_out7 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random", hclust.method = "centroid")
sel_hclust_random_out7
#> $I
#>  [1] "TMe-3164" "TMe-41"   "TMe-132"  "TMe-299"  "TMe-410"  "TMe-500" 
#>  [7] "TMe-566"  "TMe-569"  "TMe-815"  "TMe-867"  "TMe-1360" "TMe-1425"
#> [13] "TMe-1466" "TMe-1564" "TMe-1786" "TMe-1922" "TMe-2027" "TMe-2069"
#> [19] "TMe-2372" "TMe-2773" "TMe-2785" "TMe-2945" "TMe-2967" "TMe-3030"
#> [25] "TMe-3151" "TMe-3262" "TMe-3272" "TMe-3330" "TMe-3452" "TMe-3481"
#> [31] "TMe-264"  "TMe-606"  "TMe-717"  "TMe-1218" "TMe-1532" "TMe-1589"
#> [37] "TMe-2927" "TMe-2934" "TMe-2984" "TMe-2985" "TMe-2993" "TMe-3112"
#> [43] "TMe-3163" "TMe-3282" "TMe-3292" "TMe-3389" "TMe-3415" "TMe-3433"
#> [49] "TMe-3437" "TMe-3457" "TMe-3471" "TMe-3475" "TMe-3485" "TMe-3488"
#> [55] "TMe-3705" "TMe-2940" "TMe-2943" "TMe-3249" "TMe-3462" "TMe-3478"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2085" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-509"  "TMe-539" 
#>  [7] "TMe-902"  "TMe-1385" "TMe-1444" "TMe-1907" "TMe-1919" "TMe-2033"
#> [13] "TMe-2054" "TMe-2056" "TMe-2128" "TMe-2204" "TMe-2891" "TMe-2970"
#> [19] "TMe-3009" "TMe-3140" "TMe-3766" "TMe-3466" "TMe-3605" "TMe-47"  
#> [25] "TMe-171"  "TMe-478"  "TMe-1619" "TMe-2352" "TMe-2568" "TMe-2952"
#> [31] "TMe-3200" "TMe-3210" "TMe-3222" "TMe-3258" "TMe-3338" "TMe-3366"
#> [37] "TMe-3530" "TMe-3690" "TMe-3771" "TMe-3800" "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-1176" "TMe-141"  "TMe-261"  "TMe-267" 
#>  [7] "TMe-381"  "TMe-405"  "TMe-425"  "TMe-785"  "TMe-946"  "TMe-1200"
#> [13] "TMe-1593" "TMe-1868" "TMe-1897" "TMe-1939" "TMe-2897" "TMe-2968"
#> [19] "TMe-3032" "TMe-3299" "TMe-3804" "TMe-35"   "TMe-116"  "TMe-174" 
#> [25] "TMe-635"  "TMe-878"  "TMe-2405" "TMe-2756" "TMe-2926" "TMe-3143"
#> [31] "TMe-3274" "TMe-3383" "TMe-3445" "TMe-3551" "TMe-3565" "TMe-3575"
#> [37] "TMe-3715"
#> 
#> $IV
#>  [1] "TMe-3541" "TMe-154"  "TMe-274"  "TMe-388"  "TMe-427"  "TMe-184" 
#>  [7] "TMe-650"  "TMe-656"  "TMe-766"  "TMe-812"  "TMe-897"  "TMe-919" 
#> [13] "TMe-962"  "TMe-1162" "TMe-1348" "TMe-1368" "TMe-1376" "TMe-1406"
#> [19] "TMe-1434" "TMe-1525" "TMe-1806" "TMe-2025" "TMe-2043" "TMe-2059"
#> [25] "TMe-2788" "TMe-3278" "TMe-2956" "TMe-2971" "TMe-2979" "TMe-3025"
#> [31] "TMe-3054" "TMe-3055" "TMe-3089" "TMe-3255" "TMe-3256" "TMe-3273"
#> [37] "TMe-3428" "TMe-3545" "TMe-3576" "TMe-3752" "TMe-3040" "TMe-180" 
#> [43] "TMe-186"  "TMe-205"  "TMe-209"  "TMe-226"  "TMe-259"  "TMe-280" 
#> [49] "TMe-372"  "TMe-398"  "TMe-421"  "TMe-608"  "TMe-619"  "TMe-663" 
#> [55] "TMe-698"  "TMe-737"  "TMe-1020" "TMe-1129" "TMe-1511" "TMe-1988"
#> [61] "TMe-2567" "TMe-2922" "TMe-3167" "TMe-3189" "TMe-3225" "TMe-3297"
#> [67] "TMe-3386" "TMe-3417" "TMe-3422" "TMe-3527" "TMe-3701" "TMe-3707"
#> [73] "TMe-3727" "TMe-3730" "TMe-2960" "TMe-3053" "TMe-3316" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-1551" "TMe-336"  "TMe-362"  "TMe-378"  "TMe-399" 
#>  [7] "TMe-418"  "TMe-419"  "TMe-423"  "TMe-547"  "TMe-603"  "TMe-782" 
#> [13] "TMe-861"  "TMe-972"  "TMe-1099" "TMe-1204" "TMe-1215" "TMe-1366"
#> [19] "TMe-1375" "TMe-1381" "TMe-1427" "TMe-1458" "TMe-1482" "TMe-1629"
#> [25] "TMe-1924" "TMe-1963" "TMe-2041" "TMe-2050" "TMe-2060" "TMe-2121"
#> [31] "TMe-2290" "TMe-2413" "TMe-3628" "TMe-2531" "TMe-2589" "TMe-2688"
#> [37] "TMe-2790" "TMe-2853" "TMe-2863" "TMe-2973" "TMe-3185" "TMe-189" 
#> [43] "TMe-208"  "TMe-284"  "TMe-481"  "TMe-536"  "TMe-616"  "TMe-655" 
#> [49] "TMe-700"  "TMe-707"  "TMe-712"  "TMe-723"  "TMe-748"  "TMe-759" 
#> [55] "TMe-798"  "TMe-803"  "TMe-816"  "TMe-866"  "TMe-877"  "TMe-975" 
#> [61] "TMe-1042" "TMe-1196" "TMe-1199" "TMe-1305" "TMe-1311" "TMe-1332"
#> [67] "TMe-2003" "TMe-2139" "TMe-2355" "TMe-2905" "TMe-2953" "TMe-2959"
#> [73] "TMe-3034" "TMe-3311" "TMe-3356" "TMe-3408" "TMe-3736" "TMe-3659"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-130"  "TMe-1816" "TMe-269"  "TMe-361"  "TMe-442"  "TMe-465" 
#>  [7] "TMe-725"  "TMe-751"  "TMe-1383" "TMe-1403" "TMe-1416" "TMe-1428"
#> [13] "TMe-1445" "TMe-1506" "TMe-1875" "TMe-2035" "TMe-2064" "TMe-2067"
#> [19] "TMe-2196" "TMe-2534" "TMe-2791" "TMe-236"  "TMe-310"  "TMe-693" 
#> [25] "TMe-1053" "TMe-1124" "TMe-1174" "TMe-1182" "TMe-1232" "TMe-1769"
#> [31] "TMe-2510" "TMe-2543" "TMe-3095" "TMe-3387" "TMe-2983" "TMe-3116"
#> [37] "TMe-2551"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out7, use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "centroid")


## Ward's D2
sel_hclust_random_out8 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.random", hclust.method = "ward.D2")
sel_hclust_random_out8
#> $I
#>  [1] "TMe-642"  "TMe-41"   "TMe-117"  "TMe-299"  "TMe-438"  "TMe-2131"
#>  [7] "TMe-606"  "TMe-1886" "TMe-1344" "TMe-2949" "TMe-2027" "TMe-501" 
#> [13] "TMe-566"  "TMe-1448" "TMe-3480" "TMe-736"  "TMe-756"  "TMe-1533"
#> [19] "TMe-1266" "TMe-3132" "TMe-3264" "TMe-1226" "TMe-1333" "TMe-882" 
#> [25] "TMe-888"  "TMe-910"  "TMe-1777" "TMe-2372" "TMe-2031" "TMe-306" 
#> [31] "TMe-3719" "TMe-1973" "TMe-3001" "TMe-2604" "TMe-2773" "TMe-3074"
#> [37] "TMe-3211" "TMe-2975" "TMe-3548" "TMe-3478" "TMe-3481" "TMe-3726"
#> [43] "TMe-3002" "TMe-2985" "TMe-2916" "TMe-3518" "TMe-3003" "TMe-3577"
#> [49] "TMe-3641" "TMe-686"  "TMe-2984" "TMe-3501" "TMe-2927" "TMe-2911"
#> [55] "TMe-3130" "TMe-3104" "TMe-3475" "TMe-3252" "TMe-3429" "TMe-3724"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3150" "TMe-3052" "TMe-3066" "TMe-2998" "TMe-3009" "TMe-1732"
#>  [7] "TMe-3239" "TMe-902"  "TMe-1242" "TMe-2705" "TMe-2085" "TMe-455" 
#> [13] "TMe-1444" "TMe-930"  "TMe-539"  "TMe-196"  "TMe-2862" "TMe-3021"
#> [19] "TMe-660"  "TMe-1754" "TMe-1907" "TMe-2123" "TMe-2258" "TMe-2428"
#> [25] "TMe-171"  "TMe-2242" "TMe-2814" "TMe-2021" "TMe-2568" "TMe-2033"
#> [31] "TMe-1459" "TMe-2757" "TMe-74"   "TMe-2915" "TMe-3531" "TMe-3671"
#> [37] "TMe-47"   "TMe-2352" "TMe-3605" "TMe-3338" "TMe-3530"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-3029" "TMe-2977" "TMe-304"  "TMe-3544" "TMe-3048"
#>  [7] "TMe-234"  "TMe-1443" "TMe-203"  "TMe-2167" "TMe-3407" "TMe-267" 
#> [13] "TMe-420"  "TMe-1868" "TMe-773"  "TMe-2405" "TMe-946"  "TMe-1176"
#> [19] "TMe-2270" "TMe-1262" "TMe-3234" "TMe-3592" "TMe-3287" "TMe-3721"
#> [25] "TMe-2811" "TMe-3088" "TMe-2205" "TMe-3565" "TMe-3147" "TMe-3207"
#> [31] "TMe-3032" "TMe-3148" "TMe-3383" "TMe-3715" "TMe-23"   "TMe-3317"
#> [37] "TMe-3575"
#> 
#> $IV
#>  [1] "TMe-3073" "TMe-3238" "TMe-218"  "TMe-3000" "TMe-3055" "TMe-72"  
#>  [7] "TMe-108"  "TMe-150"  "TMe-2839" "TMe-460"  "TMe-318"  "TMe-266" 
#> [13] "TMe-215"  "TMe-619"  "TMe-352"  "TMe-3760" "TMe-383"  "TMe-386" 
#> [19] "TMe-387"  "TMe-2841" "TMe-3250" "TMe-3763" "TMe-1765" "TMe-556" 
#> [25] "TMe-1579" "TMe-3390" "TMe-526"  "TMe-540"  "TMe-184"  "TMe-575" 
#> [31] "TMe-641"  "TMe-1479" "TMe-675"  "TMe-221"  "TMe-3019" "TMe-834" 
#> [37] "TMe-887"  "TMe-3422" "TMe-2247" "TMe-581"  "TMe-1123" "TMe-3316"
#> [43] "TMe-1923" "TMe-761"  "TMe-1297" "TMe-1342" "TMe-3593" "TMe-1350"
#> [49] "TMe-190"  "TMe-3585" "TMe-1374" "TMe-1434" "TMe-1511" "TMe-1575"
#> [55] "TMe-2807" "TMe-3242" "TMe-1988" "TMe-2928" "TMe-2240" "TMe-2809"
#> [61] "TMe-210"  "TMe-3191" "TMe-3545" "TMe-3231" "TMe-3004" "TMe-3752"
#> [67] "TMe-3089" "TMe-3291" "TMe-3196" "TMe-3619" "TMe-3701" "TMe-3253"
#> [73] "TMe-3040" "TMe-3312" "TMe-43"   "TMe-153"  "TMe-259"  "TMe-3097"
#> [79] "TMe-1278" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-98"   "TMe-2821" "TMe-1534" "TMe-1290" "TMe-527" 
#>  [7] "TMe-294"  "TMe-932"  "TMe-1183" "TMe-2953" "TMe-1418" "TMe-197" 
#> [13] "TMe-360"  "TMe-777"  "TMe-1414" "TMe-1098" "TMe-376"  "TMe-1099"
#> [19] "TMe-1004" "TMe-391"  "TMe-390"  "TMe-1011" "TMe-417"  "TMe-288" 
#> [25] "TMe-723"  "TMe-2121" "TMe-1257" "TMe-430"  "TMe-1609" "TMe-924" 
#> [31] "TMe-481"  "TMe-1373" "TMe-2845" "TMe-1522" "TMe-600"  "TMe-603" 
#> [37] "TMe-618"  "TMe-667"  "TMe-730"  "TMe-819"  "TMe-861"  "TMe-1390"
#> [43] "TMe-1877" "TMe-2339" "TMe-2835" "TMe-813"  "TMe-1026" "TMe-1042"
#> [49] "TMe-1523" "TMe-803"  "TMe-647"  "TMe-2933" "TMe-1366" "TMe-1427"
#> [55] "TMe-343"  "TMe-997"  "TMe-1643" "TMe-1963" "TMe-181"  "TMe-742" 
#> [61] "TMe-1220" "TMe-1391" "TMe-1934" "TMe-2015" "TMe-2050" "TMe-2139"
#> [67] "TMe-759"  "TMe-2151" "TMe-3034" "TMe-2761" "TMe-2435" "TMe-3185"
#> [73] "TMe-1293" "TMe-344"  "TMe-687"  "TMe-748"  "TMe-1733" "TMe-2905"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-1182" "TMe-1985" "TMe-2398" "TMe-1007" "TMe-361"  "TMe-1383"
#>  [7] "TMe-1902" "TMe-1033" "TMe-548"  "TMe-1509" "TMe-598"  "TMe-3177"
#> [13] "TMe-751"  "TMe-2791" "TMe-1074" "TMe-846"  "TMe-1238" "TMe-1676"
#> [19] "TMe-1460" "TMe-1261" "TMe-1178" "TMe-1483" "TMe-1124" "TMe-1428"
#> [25] "TMe-693"  "TMe-963"  "TMe-1506" "TMe-3368" "TMe-1992" "TMe-1832"
#> [31] "TMe-1875" "TMe-2035" "TMe-1744" "TMe-2196" "TMe-188"  "TMe-1769"
#> [37] "TMe-3549"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_random_out8, use.names = FALSE)) +
  labs(title = "hclust.random", subtitle = "ward.D2")


# Cluster-Based Sampling via Hierarchical Clustering with Medoid Selection

## UPGMA
sel_hclust_medoid_out1 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid", hclust.method = "average")
sel_hclust_medoid_out1
#> $I
#>  [1] "TMe-1981" "TMe-41"   "TMe-2810" "TMe-1717" "TMe-132"  "TMe-3424"
#>  [7] "TMe-299"  "TMe-486"  "TMe-1564" "TMe-1823" "TMe-469"  "TMe-489" 
#> [13] "TMe-500"  "TMe-566"  "TMe-1448" "TMe-642"  "TMe-3252" "TMe-306" 
#> [19] "TMe-3132" "TMe-815"  "TMe-867"  "TMe-610"  "TMe-882"  "TMe-1086"
#> [25] "TMe-3441" "TMe-1964" "TMe-2253" "TMe-3490" "TMe-1777" "TMe-1940"
#> [31] "TMe-2010" "TMe-2604" "TMe-3477" "TMe-2967" "TMe-3729" "TMe-2964"
#> [37] "TMe-3235" "TMe-3151" "TMe-3169" "TMe-3419" "TMe-3272" "TMe-3330"
#> [43] "TMe-2913" "TMe-2917" "TMe-3685" "TMe-2462" "TMe-3501" "TMe-2927"
#> [49] "TMe-2984" "TMe-2996" "TMe-3499" "TMe-3026" "TMe-3163" "TMe-3292"
#> [55] "TMe-3389" "TMe-3415" "TMe-3437" "TMe-3471" "TMe-3694" "TMe-3478"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-40"   "TMe-2611" "TMe-842"  "TMe-1250" "TMe-3009" "TMe-289" 
#>  [7] "TMe-3617" "TMe-369"  "TMe-2123" "TMe-1969" "TMe-3203" "TMe-2329"
#> [13] "TMe-2951" "TMe-2824" "TMe-902"  "TMe-1005" "TMe-2414" "TMe-1385"
#> [19] "TMe-2021" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2056" "TMe-3498"
#> [25] "TMe-2851" "TMe-3455" "TMe-3642" "TMe-3101" "TMe-3766" "TMe-3466"
#> [31] "TMe-3605" "TMe-176"  "TMe-1668" "TMe-2952" "TMe-3200" "TMe-3210"
#> [37] "TMe-3222" "TMe-3237" "TMe-3338" "TMe-3406" "TMe-3690"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-1443" "TMe-2977" "TMe-93"   "TMe-3048" "TMe-123" 
#>  [7] "TMe-141"  "TMe-187"  "TMe-2161" "TMe-206"  "TMe-1138" "TMe-261" 
#> [13] "TMe-267"  "TMe-405"  "TMe-420"  "TMe-434"  "TMe-2394" "TMe-1787"
#> [19] "TMe-1176" "TMe-3324" "TMe-2356" "TMe-1939" "TMe-3100" "TMe-2205"
#> [25] "TMe-2897" "TMe-3663" "TMe-3592" "TMe-3094" "TMe-3128" "TMe-3750"
#> [31] "TMe-3299" "TMe-3638" "TMe-3721" "TMe-3598" "TMe-35"   "TMe-2926"
#> [37] "TMe-3143"
#> 
#> $IV
#>  [1] "TMe-12"   "TMe-3068" "TMe-3780" "TMe-1374" "TMe-107"  "TMe-634" 
#>  [7] "TMe-1364" "TMe-2839" "TMe-318"  "TMe-415"  "TMe-802"  "TMe-314" 
#> [13] "TMe-2954" "TMe-82"   "TMe-218"  "TMe-3537" "TMe-386"  "TMe-552" 
#> [19] "TMe-444"  "TMe-394"  "TMe-396"  "TMe-1726" "TMe-427"  "TMe-3168"
#> [25] "TMe-2172" "TMe-513"  "TMe-2039" "TMe-550"  "TMe-3542" "TMe-1470"
#> [31] "TMe-1139" "TMe-1489" "TMe-812"  "TMe-897"  "TMe-919"  "TMe-1336"
#> [37] "TMe-1377" "TMe-2567" "TMe-2399" "TMe-1350" "TMe-190"  "TMe-3409"
#> [43] "TMe-1434" "TMe-153"  "TMe-1150" "TMe-2928" "TMe-2043" "TMe-2240"
#> [49] "TMe-78"   "TMe-2809" "TMe-2833" "TMe-2971" "TMe-3538" "TMe-2989"
#> [55] "TMe-3025" "TMe-3054" "TMe-3055" "TMe-3594" "TMe-3255" "TMe-3273"
#> [61] "TMe-3253" "TMe-3283" "TMe-3428" "TMe-3545" "TMe-3758" "TMe-3581"
#> [67] "TMe-3494" "TMe-2758" "TMe-714"  "TMe-3602" "TMe-698"  "TMe-737" 
#> [73] "TMe-1328" "TMe-1020" "TMe-3109" "TMe-1511" "TMe-3730" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-99"   "TMe-1127" "TMe-2769" "TMe-168"  "TMe-447" 
#>  [7] "TMe-294"  "TMe-629"  "TMe-1160" "TMe-334"  "TMe-951"  "TMe-3332"
#> [13] "TMe-1446" "TMe-2518" "TMe-1234" "TMe-1011" "TMe-1067" "TMe-378" 
#> [19] "TMe-385"  "TMe-390"  "TMe-800"  "TMe-418"  "TMe-1455" "TMe-2032"
#> [25] "TMe-982"  "TMe-932"  "TMe-476"  "TMe-645"  "TMe-884"  "TMe-603" 
#> [31] "TMe-730"  "TMe-819"  "TMe-828"  "TMe-861"  "TMe-1411" "TMe-1879"
#> [37] "TMe-1295" "TMe-813"  "TMe-1099" "TMe-1357" "TMe-1307" "TMe-1440"
#> [43] "TMe-1482" "TMe-1500" "TMe-1629" "TMe-2542" "TMe-2753" "TMe-1871"
#> [49] "TMe-2192" "TMe-1733" "TMe-2041" "TMe-2050" "TMe-2121" "TMe-2290"
#> [55] "TMe-1487" "TMe-2413" "TMe-2435" "TMe-2530" "TMe-2589" "TMe-2790"
#> [61] "TMe-2853" "TMe-3185" "TMe-189"  "TMe-771"  "TMe-284"  "TMe-481" 
#> [67] "TMe-536"  "TMe-712"  "TMe-723"  "TMe-750"  "TMe-1042" "TMe-1196"
#> [73] "TMe-2355" "TMe-2905" "TMe-2961" "TMe-2959" "TMe-3356" "TMe-3736"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-1985" "TMe-1017" "TMe-1477" "TMe-1512" "TMe-442" 
#>  [7] "TMe-465"  "TMe-626"  "TMe-1902" "TMe-548"  "TMe-1110" "TMe-791" 
#> [13] "TMe-725"  "TMe-281"  "TMe-876"  "TMe-1079" "TMe-2963" "TMe-1832"
#> [19] "TMe-2389" "TMe-1383" "TMe-1413" "TMe-1416" "TMe-1661" "TMe-1445"
#> [25] "TMe-1518" "TMe-1775" "TMe-1875" "TMe-2035" "TMe-2196" "TMe-2791"
#> [31] "TMe-188"  "TMe-1053" "TMe-1174" "TMe-1646" "TMe-1769" "TMe-2510"
#> [37] "TMe-3095"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out1, use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "average")


## Single-linkage
sel_hclust_medoid_out2 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid", hclust.method = "single")
sel_hclust_medoid_out2
#> $I
#>  [1] "TMe-1696" "TMe-41"   "TMe-132"  "TMe-469"  "TMe-500"  "TMe-566" 
#>  [7] "TMe-569"  "TMe-642"  "TMe-815"  "TMe-867"  "TMe-1096" "TMe-1425"
#> [13] "TMe-1466" "TMe-1564" "TMe-1716" "TMe-2027" "TMe-2040" "TMe-2069"
#> [19] "TMe-2372" "TMe-2773" "TMe-2785" "TMe-2967" "TMe-3151" "TMe-3169"
#> [25] "TMe-3262" "TMe-3272" "TMe-3345" "TMe-3353" "TMe-3452" "TMe-3633"
#> [31] "TMe-606"  "TMe-715"  "TMe-717"  "TMe-1218" "TMe-2462" "TMe-2912"
#> [37] "TMe-2934" "TMe-2937" "TMe-2984" "TMe-2985" "TMe-2993" "TMe-3112"
#> [43] "TMe-3163" "TMe-3229" "TMe-3282" "TMe-3292" "TMe-3389" "TMe-3415"
#> [49] "TMe-3433" "TMe-3437" "TMe-3471" "TMe-3475" "TMe-3488" "TMe-3514"
#> [55] "TMe-3694" "TMe-2940" "TMe-2943" "TMe-3462" "TMe-3463" "TMe-3478"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-1969" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-404"  "TMe-768" 
#>  [7] "TMe-902"  "TMe-1242" "TMe-1385" "TMe-1444" "TMe-1831" "TMe-1907"
#> [13] "TMe-1919" "TMe-2033" "TMe-2054" "TMe-2056" "TMe-2204" "TMe-2851"
#> [19] "TMe-2866" "TMe-2891" "TMe-3009" "TMe-3140" "TMe-3141" "TMe-3466"
#> [25] "TMe-3605" "TMe-47"   "TMe-176"  "TMe-431"  "TMe-478"  "TMe-648" 
#> [31] "TMe-1310" "TMe-2352" "TMe-2568" "TMe-2952" "TMe-3200" "TMe-3210"
#> [37] "TMe-3222" "TMe-3338" "TMe-3690" "TMe-3771" "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-2154" "TMe-104"  "TMe-261"  "TMe-267" 
#>  [7] "TMe-381"  "TMe-425"  "TMe-785"  "TMe-1200" "TMe-1593" "TMe-1897"
#> [13] "TMe-2270" "TMe-2374" "TMe-2823" "TMe-2897" "TMe-3299" "TMe-3804"
#> [19] "TMe-35"   "TMe-116"  "TMe-174"  "TMe-187"  "TMe-234"  "TMe-635" 
#> [25] "TMe-2756" "TMe-3598" "TMe-3143" "TMe-3274" "TMe-3383" "TMe-3397"
#> [31] "TMe-3445" "TMe-3551" "TMe-3565" "TMe-3575" "TMe-3596" "TMe-3620"
#> [37] "TMe-3631"
#> 
#> $IV
#>  [1] "TMe-2340" "TMe-154"  "TMe-314"  "TMe-427"  "TMe-650"  "TMe-766" 
#>  [7] "TMe-812"  "TMe-897"  "TMe-919"  "TMe-956"  "TMe-962"  "TMe-967" 
#> [13] "TMe-1342" "TMe-1350" "TMe-1368" "TMe-1369" "TMe-1376" "TMe-1406"
#> [19] "TMe-1434" "TMe-1525" "TMe-1806" "TMe-2020" "TMe-2025" "TMe-2043"
#> [25] "TMe-2809" "TMe-2839" "TMe-2956" "TMe-2971" "TMe-3025" "TMe-3054"
#> [31] "TMe-3055" "TMe-3089" "TMe-3255" "TMe-3256" "TMe-3273" "TMe-3275"
#> [37] "TMe-3283" "TMe-3545" "TMe-5"    "TMe-25"   "TMe-180"  "TMe-205" 
#> [43] "TMe-209"  "TMe-225"  "TMe-226"  "TMe-259"  "TMe-280"  "TMe-372" 
#> [49] "TMe-608"  "TMe-619"  "TMe-663"  "TMe-698"  "TMe-733"  "TMe-735" 
#> [55] "TMe-737"  "TMe-1010" "TMe-1020" "TMe-1129" "TMe-1148" "TMe-1181"
#> [61] "TMe-1511" "TMe-1988" "TMe-2567" "TMe-2922" "TMe-3167" "TMe-3189"
#> [67] "TMe-3225" "TMe-3228" "TMe-3297" "TMe-3386" "TMe-3422" "TMe-3527"
#> [73] "TMe-3701" "TMe-3707" "TMe-3730" "TMe-3053" "TMe-3316" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-332"  "TMe-347"  "TMe-376"  "TMe-378"  "TMe-399" 
#>  [7] "TMe-418"  "TMe-419"  "TMe-423"  "TMe-782"  "TMe-819"  "TMe-861" 
#> [13] "TMe-972"  "TMe-1099" "TMe-1204" "TMe-1357" "TMe-1366" "TMe-1381"
#> [19] "TMe-1458" "TMe-1482" "TMe-1488" "TMe-1730" "TMe-1873" "TMe-1924"
#> [25] "TMe-1963" "TMe-2041" "TMe-2050" "TMe-2060" "TMe-2121" "TMe-2192"
#> [31] "TMe-2290" "TMe-2413" "TMe-2531" "TMe-2589" "TMe-2612" "TMe-2688"
#> [37] "TMe-2750" "TMe-2790" "TMe-2853" "TMe-2863" "TMe-2973" "TMe-3185"
#> [43] "TMe-189"  "TMe-194"  "TMe-208"  "TMe-277"  "TMe-284"  "TMe-341" 
#> [49] "TMe-481"  "TMe-536"  "TMe-616"  "TMe-700"  "TMe-712"  "TMe-723" 
#> [55] "TMe-832"  "TMe-866"  "TMe-877"  "TMe-1042" "TMe-1196" "TMe-1199"
#> [61] "TMe-1300" "TMe-1305" "TMe-1307" "TMe-1311" "TMe-1391" "TMe-1453"
#> [67] "TMe-2003" "TMe-2139" "TMe-2355" "TMe-2905" "TMe-2953" "TMe-2959"
#> [73] "TMe-3034" "TMe-3311" "TMe-3356" "TMe-3363" "TMe-3408" "TMe-3736"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-130"  "TMe-2195" "TMe-269"  "TMe-1550" "TMe-361"  "TMe-442" 
#>  [7] "TMe-465"  "TMe-2308" "TMe-1383" "TMe-1416" "TMe-1428" "TMe-1518"
#> [13] "TMe-1775" "TMe-1875" "TMe-2035" "TMe-2064" "TMe-2067" "TMe-2196"
#> [19] "TMe-2534" "TMe-690"  "TMe-693"  "TMe-1053" "TMe-1124" "TMe-1174"
#> [25] "TMe-1182" "TMe-1232" "TMe-1239" "TMe-1646" "TMe-1723" "TMe-1769"
#> [31] "TMe-2510" "TMe-2963" "TMe-3095" "TMe-3387" "TMe-2983" "TMe-3116"
#> [37] "TMe-3549"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out2, use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "single")


## Complete-linkage
sel_hclust_medoid_out3 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid", hclust.method = "complete")
sel_hclust_medoid_out3
#> $I
#>  [1] "TMe-642"  "TMe-1621" "TMe-33"   "TMe-2383" "TMe-3480" "TMe-3424"
#>  [7] "TMe-1981" "TMe-2131" "TMe-1564" "TMe-1930" "TMe-882"  "TMe-500" 
#> [13] "TMe-1696" "TMe-1190" "TMe-2083" "TMe-3170" "TMe-3633" "TMe-3252"
#> [19] "TMe-306"  "TMe-28"   "TMe-815"  "TMe-3266" "TMe-990"  "TMe-610" 
#> [25] "TMe-3127" "TMe-1086" "TMe-1341" "TMe-910"  "TMe-1964" "TMe-20"  
#> [31] "TMe-1955" "TMe-1940" "TMe-1532" "TMe-3392" "TMe-2810" "TMe-2917"
#> [37] "TMe-2936" "TMe-3065" "TMe-3262" "TMe-2913" "TMe-3416" "TMe-3211"
#> [43] "TMe-3296" "TMe-3521" "TMe-3151" "TMe-3220" "TMe-2985" "TMe-3330"
#> [49] "TMe-2909" "TMe-3236" "TMe-3685" "TMe-248"  "TMe-3501" "TMe-2927"
#> [55] "TMe-2996" "TMe-3419" "TMe-3398" "TMe-3425" "TMe-3694" "TMe-2943"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-2"    "TMe-2611" "TMe-3571" "TMe-85"   "TMe-3009" "TMe-2414"
#>  [7] "TMe-377"  "TMe-2715" "TMe-2765" "TMe-433"  "TMe-1833" "TMe-455" 
#> [13] "TMe-477"  "TMe-589"  "TMe-2951" "TMe-1484" "TMe-1505" "TMe-2824"
#> [19] "TMe-3365" "TMe-1864" "TMe-674"  "TMe-2123" "TMe-1250" "TMe-176" 
#> [25] "TMe-3146" "TMe-2021" "TMe-1271" "TMe-2033" "TMe-2056" "TMe-2995"
#> [31] "TMe-3498" "TMe-3455" "TMe-3642" "TMe-3093" "TMe-3101" "TMe-3766"
#> [37] "TMe-74"   "TMe-706"  "TMe-2952" "TMe-3188" "TMe-3284"
#> 
#> $III
#>  [1] "TMe-155"  "TMe-1443" "TMe-14"   "TMe-93"   "TMe-995"  "TMe-123" 
#>  [7] "TMe-1262" "TMe-161"  "TMe-1797" "TMe-3234" "TMe-267"  "TMe-420" 
#> [13] "TMe-434"  "TMe-1738" "TMe-2394" "TMe-946"  "TMe-3321" "TMe-1804"
#> [19] "TMe-1819" "TMe-1965" "TMe-1421" "TMe-2169" "TMe-3643" "TMe-3324"
#> [25] "TMe-3207" "TMe-3750" "TMe-3043" "TMe-617"  "TMe-3128" "TMe-3299"
#> [31] "TMe-3638" "TMe-3804" "TMe-23"   "TMe-187"  "TMe-3143" "TMe-3445"
#> [37] "TMe-3565"
#> 
#> $IV
#>  [1] "TMe-3073" "TMe-3045" "TMe-57"   "TMe-3084" "TMe-103"  "TMe-108" 
#>  [7] "TMe-1364" "TMe-154"  "TMe-170"  "TMe-534"  "TMe-834"  "TMe-266" 
#> [13] "TMe-65"   "TMe-314"  "TMe-318"  "TMe-2954" "TMe-552"  "TMe-82"  
#> [19] "TMe-218"  "TMe-1167" "TMe-386"  "TMe-444"  "TMe-1148" "TMe-396" 
#> [25] "TMe-3144" "TMe-1336" "TMe-519"  "TMe-3168" "TMe-2039" "TMe-3435"
#> [31] "TMe-87"   "TMe-550"  "TMe-3212" "TMe-1402" "TMe-1700" "TMe-1139"
#> [37] "TMe-1430" "TMe-3161" "TMe-3019" "TMe-3541" "TMe-802"  "TMe-1923"
#> [43] "TMe-1891" "TMe-699"  "TMe-2980" "TMe-3585" "TMe-885"  "TMe-1027"
#> [49] "TMe-1726" "TMe-1150" "TMe-2210" "TMe-1903" "TMe-3231" "TMe-2928"
#> [55] "TMe-708"  "TMe-2833" "TMe-3357" "TMe-2809" "TMe-2971" "TMe-3494"
#> [61] "TMe-3068" "TMe-3198" "TMe-3025" "TMe-3053" "TMe-3594" "TMe-3255"
#> [67] "TMe-3253" "TMe-3283" "TMe-3781" "TMe-27"   "TMe-2775" "TMe-2947"
#> [73] "TMe-2758" "TMe-1419" "TMe-698"  "TMe-1020" "TMe-3097" "TMe-3297"
#> [79] "TMe-3417" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-208"  "TMe-969"  "TMe-1127" "TMe-167"  "TMe-212"  "TMe-868" 
#>  [7] "TMe-294"  "TMe-1101" "TMe-307"  "TMe-2441" "TMe-1399" "TMe-2213"
#> [13] "TMe-3332" "TMe-1358" "TMe-1680" "TMe-948"  "TMe-982"  "TMe-1099"
#> [19] "TMe-629"  "TMe-385"  "TMe-390"  "TMe-1011" "TMe-402"  "TMe-424" 
#> [25] "TMe-645"  "TMe-474"  "TMe-1455" "TMe-895"  "TMe-1834" "TMe-870" 
#> [31] "TMe-574"  "TMe-1299" "TMe-543"  "TMe-1234" "TMe-623"  "TMe-730" 
#> [37] "TMe-777"  "TMe-1227" "TMe-2138" "TMe-861"  "TMe-2835" "TMe-1523"
#> [43] "TMe-1357" "TMe-1307" "TMe-1901" "TMe-137"  "TMe-2009" "TMe-1963"
#> [49] "TMe-1500" "TMe-997"  "TMe-1629" "TMe-2753" "TMe-1871" "TMe-1391"
#> [55] "TMe-1269" "TMe-2050" "TMe-2121" "TMe-2151" "TMe-2290" "TMe-1487"
#> [61] "TMe-2435" "TMe-2530" "TMe-2589" "TMe-2853" "TMe-3185" "TMe-1188"
#> [67] "TMe-189"  "TMe-764"  "TMe-284"  "TMe-685"  "TMe-341"  "TMe-397" 
#> [73] "TMe-481"  "TMe-536"  "TMe-723"  "TMe-1733" "TMe-1785" "TMe-2961"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-1985" "TMe-1017" "TMe-1477" "TMe-1512" "TMe-442" 
#>  [7] "TMe-465"  "TMe-626"  "TMe-1902" "TMe-531"  "TMe-1110" "TMe-2402"
#> [13] "TMe-281"  "TMe-705"  "TMe-1164" "TMe-838"  "TMe-876"  "TMe-1079"
#> [19] "TMe-1460" "TMe-2963" "TMe-1832" "TMe-625"  "TMe-3095" "TMe-1413"
#> [25] "TMe-1416" "TMe-1661" "TMe-1445" "TMe-2067" "TMe-1518" "TMe-2389"
#> [31] "TMe-1995" "TMe-2035" "TMe-1608" "TMe-188"  "TMe-1053" "TMe-1174"
#> [37] "TMe-1646"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out3, use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "complete")


## Ward's D
sel_hclust_medoid_out4 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid", hclust.method = "ward.D")
sel_hclust_medoid_out4
#> $I
#>  [1] "TMe-642"  "TMe-2810" "TMe-3149" "TMe-300"  "TMe-3424" "TMe-2453"
#>  [7] "TMe-1564" "TMe-489"  "TMe-3202" "TMe-306"  "TMe-677"  "TMe-3441"
#> [13] "TMe-1448" "TMe-1981" "TMe-1341" "TMe-756"  "TMe-2253" "TMe-3521"
#> [19] "TMe-3132" "TMe-3330" "TMe-1823" "TMe-610"  "TMe-882"  "TMe-1170"
#> [25] "TMe-1091" "TMe-1086" "TMe-2131" "TMe-1964" "TMe-486"  "TMe-910" 
#> [31] "TMe-1621" "TMe-1869" "TMe-1940" "TMe-2066" "TMe-2152" "TMe-3223"
#> [37] "TMe-3477" "TMe-3490" "TMe-3027" "TMe-3262" "TMe-3296" "TMe-3433"
#> [43] "TMe-3220" "TMe-3268" "TMe-2985" "TMe-2916" "TMe-2917" "TMe-3236"
#> [49] "TMe-3685" "TMe-248"  "TMe-3698" "TMe-3501" "TMe-3425" "TMe-3493"
#> [55] "TMe-2975" "TMe-3463" "TMe-3102" "TMe-3398" "TMe-3170" "TMe-3252"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-40"   "TMe-2611" "TMe-3571" "TMe-85"   "TMe-3009" "TMe-2414"
#>  [7] "TMe-377"  "TMe-2715" "TMe-478"  "TMe-2765" "TMe-2950" "TMe-1474"
#> [13] "TMe-455"  "TMe-2329" "TMe-706"  "TMe-2951" "TMe-589"  "TMe-1484"
#> [19] "TMe-2843" "TMe-2824" "TMe-2200" "TMe-674"  "TMe-2123" "TMe-2258"
#> [25] "TMe-1698" "TMe-176"  "TMe-1833" "TMe-1505" "TMe-2021" "TMe-2033"
#> [31] "TMe-2056" "TMe-3498" "TMe-74"   "TMe-3455" "TMe-3642" "TMe-3093"
#> [37] "TMe-2995" "TMe-3540" "TMe-3766" "TMe-2952" "TMe-3188"
#> 
#> $III
#>  [1] "TMe-161"  "TMe-3029" "TMe-3644" "TMe-1836" "TMe-3302" "TMe-3048"
#>  [7] "TMe-187"  "TMe-200"  "TMe-206"  "TMe-1138" "TMe-3005" "TMe-267" 
#> [13] "TMe-420"  "TMe-434"  "TMe-853"  "TMe-2394" "TMe-3028" "TMe-1176"
#> [19] "TMe-1398" "TMe-3118" "TMe-3277" "TMe-1207" "TMe-2158" "TMe-1288"
#> [25] "TMe-3643" "TMe-3324" "TMe-3382" "TMe-3750" "TMe-3043" "TMe-3592"
#> [31] "TMe-3128" "TMe-3638" "TMe-3598" "TMe-3317" "TMe-3143" "TMe-3572"
#> [37] "TMe-3631"
#> 
#> $IV
#>  [1] "TMe-3045" "TMe-204"  "TMe-62"   "TMe-1765" "TMe-72"   "TMe-634" 
#>  [7] "TMe-1364" "TMe-2839" "TMe-207"  "TMe-318"  "TMe-270"  "TMe-594" 
#> [13] "TMe-2954" "TMe-82"   "TMe-1406" "TMe-386"  "TMe-387"  "TMe-1148"
#> [19] "TMe-396"  "TMe-1726" "TMe-3144" "TMe-3168" "TMe-1579" "TMe-513" 
#> [25] "TMe-3435" "TMe-87"   "TMe-550"  "TMe-3542" "TMe-840"  "TMe-1700"
#> [31] "TMe-1144" "TMe-3757" "TMe-3019" "TMe-834"  "TMe-891"  "TMe-18"  
#> [37] "TMe-1139" "TMe-1402" "TMe-2039" "TMe-3727" "TMe-1151" "TMe-1419"
#> [43] "TMe-1167" "TMe-802"  "TMe-1342" "TMe-2399" "TMe-1350" "TMe-190" 
#> [49] "TMe-3409" "TMe-885"  "TMe-1027" "TMe-1485" "TMe-1336" "TMe-3537"
#> [55] "TMe-3451" "TMe-1903" "TMe-2928" "TMe-209"  "TMe-2210" "TMe-2958"
#> [61] "TMe-2809" "TMe-3240" "TMe-3615" "TMe-3025" "TMe-3053" "TMe-3073"
#> [67] "TMe-3594" "TMe-3161" "TMe-3108" "TMe-3562" "TMe-3257" "TMe-3758"
#> [73] "TMe-3591" "TMe-3581" "TMe-3044" "TMe-63"   "TMe-2758" "TMe-1278"
#> [79] "TMe-3253" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-99"   "TMe-1127" "TMe-1534" "TMe-212"  "TMe-275" 
#>  [7] "TMe-294"  "TMe-305"  "TMe-1299" "TMe-1399" "TMe-2009" "TMe-3332"
#> [13] "TMe-1358" "TMe-2441" "TMe-826"  "TMe-1234" "TMe-623"  "TMe-1067"
#> [19] "TMe-1099" "TMe-939"  "TMe-385"  "TMe-1100" "TMe-1011" "TMe-800" 
#> [25] "TMe-1322" "TMe-527"  "TMe-217"  "TMe-1458" "TMe-1859" "TMe-430" 
#> [31] "TMe-429"  "TMe-456"  "TMe-211"  "TMe-476"  "TMe-487"  "TMe-870" 
#> [37] "TMe-1183" "TMe-982"  "TMe-2845" "TMe-884"  "TMe-667"  "TMe-730" 
#> [43] "TMe-750"  "TMe-2355" "TMe-1411" "TMe-997"  "TMe-1877" "TMe-2213"
#> [49] "TMe-972"  "TMe-813"  "TMe-1523" "TMe-2290" "TMe-1188" "TMe-2933"
#> [55] "TMe-1307" "TMe-1901" "TMe-1440" "TMe-1500" "TMe-1694" "TMe-2753"
#> [61] "TMe-1308" "TMe-1788" "TMe-1871" "TMe-2192" "TMe-858"  "TMe-1269"
#> [67] "TMe-2050" "TMe-850"  "TMe-2530" "TMe-1487" "TMe-2435" "TMe-3185"
#> [73] "TMe-669"  "TMe-284"  "TMe-1227" "TMe-771"  "TMe-1733" "TMe-2905"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-1985" "TMe-678"  "TMe-1477" "TMe-1512" "TMe-620" 
#>  [7] "TMe-2308" "TMe-626"  "TMe-548"  "TMe-1110" "TMe-2402" "TMe-273" 
#> [13] "TMe-705"  "TMe-1481" "TMe-1164" "TMe-846"  "TMe-876"  "TMe-1079"
#> [19] "TMe-281"  "TMe-1460" "TMe-2963" "TMe-625"  "TMe-1413" "TMe-1416"
#> [25] "TMe-1661" "TMe-1486" "TMe-1518" "TMe-2389" "TMe-1775" "TMe-1832"
#> [31] "TMe-2551" "TMe-2035" "TMe-2196" "TMe-188"  "TMe-1053" "TMe-1608"
#> [37] "TMe-3095"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out4, use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "ward.D")


## WPGMA
sel_hclust_medoid_out5 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid", hclust.method = "mcquitty")
sel_hclust_medoid_out5
#> $I
#>  [1] "TMe-1981" "TMe-2810" "TMe-3149" "TMe-132"  "TMe-3424" "TMe-299" 
#>  [7] "TMe-1890" "TMe-1564" "TMe-1823" "TMe-3202" "TMe-882"  "TMe-500" 
#> [13] "TMe-2964" "TMe-566"  "TMe-1448" "TMe-642"  "TMe-1341" "TMe-3252"
#> [19] "TMe-306"  "TMe-778"  "TMe-815"  "TMe-867"  "TMe-610"  "TMe-1086"
#> [25] "TMe-3441" "TMe-2383" "TMe-1964" "TMe-3490" "TMe-1777" "TMe-2103"
#> [31] "TMe-1940" "TMe-2604" "TMe-2773" "TMe-3132" "TMe-2967" "TMe-3296"
#> [37] "TMe-3151" "TMe-3220" "TMe-3577" "TMe-3272" "TMe-3330" "TMe-2917"
#> [43] "TMe-3685" "TMe-717"  "TMe-1532" "TMe-2462" "TMe-3104" "TMe-2913"
#> [49] "TMe-2927" "TMe-3460" "TMe-2984" "TMe-3026" "TMe-3163" "TMe-3292"
#> [55] "TMe-3389" "TMe-3415" "TMe-3437" "TMe-3471" "TMe-3694" "TMe-3493"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3540" "TMe-3533" "TMe-3571" "TMe-160"  "TMe-289"  "TMe-431" 
#>  [7] "TMe-369"  "TMe-2123" "TMe-1864" "TMe-2242" "TMe-2329" "TMe-589" 
#> [13] "TMe-2951" "TMe-2843" "TMe-2824" "TMe-1474" "TMe-902"  "TMe-1005"
#> [19] "TMe-2414" "TMe-3365" "TMe-1385" "TMe-1422" "TMe-176"  "TMe-1831"
#> [25] "TMe-2021" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2056" "TMe-3498"
#> [31] "TMe-3642" "TMe-3766" "TMe-3466" "TMe-3605" "TMe-1668" "TMe-3200"
#> [37] "TMe-3286" "TMe-3284" "TMe-3338" "TMe-3690" "TMe-3800"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-2502" "TMe-2977" "TMe-3069" "TMe-93"   "TMe-161" 
#>  [7] "TMe-3029" "TMe-3302" "TMe-187"  "TMe-3772" "TMe-1207" "TMe-1138"
#> [13] "TMe-261"  "TMe-267"  "TMe-405"  "TMe-420"  "TMe-3028" "TMe-853" 
#> [19] "TMe-785"  "TMe-1787" "TMe-1790" "TMe-2356" "TMe-3721" "TMe-1939"
#> [25] "TMe-2362" "TMe-2897" "TMe-3382" "TMe-3663" "TMe-3592" "TMe-3598"
#> [31] "TMe-3128" "TMe-3750" "TMe-3299" "TMe-3638" "TMe-23"   "TMe-2926"
#> [37] "TMe-3143"
#> 
#> $IV
#>  [1] "TMe-513"  "TMe-3045" "TMe-3068" "TMe-3562" "TMe-1652" "TMe-72"  
#>  [7] "TMe-634"  "TMe-1364" "TMe-2839" "TMe-318"  "TMe-270"  "TMe-314" 
#> [13] "TMe-2954" "TMe-82"   "TMe-840"  "TMe-2172" "TMe-2390" "TMe-552" 
#> [19] "TMe-444"  "TMe-1148" "TMe-3231" "TMe-1726" "TMe-427"  "TMe-1903"
#> [25] "TMe-2039" "TMe-540"  "TMe-550"  "TMe-1700" "TMe-699"  "TMe-812" 
#> [31] "TMe-897"  "TMe-919"  "TMe-235"  "TMe-1278" "TMe-1377" "TMe-3752"
#> [37] "TMe-802"  "TMe-1342" "TMe-2980" "TMe-1350" "TMe-190"  "TMe-3409"
#> [43] "TMe-885"  "TMe-1434" "TMe-153"  "TMe-1489" "TMe-280"  "TMe-1150"
#> [49] "TMe-2947" "TMe-1916" "TMe-2928" "TMe-2043" "TMe-2240" "TMe-78"  
#> [55] "TMe-2809" "TMe-2833" "TMe-3538" "TMe-2989" "TMe-3025" "TMe-3054"
#> [61] "TMe-3055" "TMe-3594" "TMe-3535" "TMe-3255" "TMe-3273" "TMe-3253"
#> [67] "TMe-3283" "TMe-3581" "TMe-2758" "TMe-3602" "TMe-698"  "TMe-737" 
#> [73] "TMe-1328" "TMe-1020" "TMe-1597" "TMe-3189" "TMe-3730" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-99"   "TMe-1127" "TMe-2769" "TMe-168"  "TMe-1405"
#>  [7] "TMe-197"  "TMe-432"  "TMe-359"  "TMe-1160" "TMe-334"  "TMe-951" 
#> [13] "TMe-3332" "TMe-362"  "TMe-474"  "TMe-1011" "TMe-378"  "TMe-629" 
#> [19] "TMe-385"  "TMe-390"  "TMe-800"  "TMe-418"  "TMe-1455" "TMe-1859"
#> [25] "TMe-982"  "TMe-475"  "TMe-476"  "TMe-1199" "TMe-2400" "TMe-645" 
#> [31] "TMe-1234" "TMe-603"  "TMe-730"  "TMe-2530" "TMe-669"  "TMe-861" 
#> [37] "TMe-1411" "TMe-926"  "TMe-1227" "TMe-972"  "TMe-1099" "TMe-1192"
#> [43] "TMe-1295" "TMe-1357" "TMe-1307" "TMe-1901" "TMe-1482" "TMe-1500"
#> [49] "TMe-1629" "TMe-2326" "TMe-2753" "TMe-1871" "TMe-2032" "TMe-1269"
#> [55] "TMe-2050" "TMe-2121" "TMe-2290" "TMe-1487" "TMe-2413" "TMe-2435"
#> [61] "TMe-2531" "TMe-2589" "TMe-2853" "TMe-2973" "TMe-3185" "TMe-771" 
#> [67] "TMe-284"  "TMe-341"  "TMe-536"  "TMe-712"  "TMe-723"  "TMe-1042"
#> [73] "TMe-1785" "TMe-2905" "TMe-2961" "TMe-2959" "TMe-3356" "TMe-3736"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-1985" "TMe-1017" "TMe-626"  "TMe-1614" "TMe-442" 
#>  [7] "TMe-465"  "TMe-2308" "TMe-548"  "TMe-1110" "TMe-2402" "TMe-281" 
#> [13] "TMe-725"  "TMe-1477" "TMe-1164" "TMe-876"  "TMe-1079" "TMe-1945"
#> [19] "TMe-2963" "TMe-625"  "TMe-620"  "TMe-1413" "TMe-1416" "TMe-1445"
#> [25] "TMe-1518" "TMe-2389" "TMe-1832" "TMe-2551" "TMe-2035" "TMe-1608"
#> [31] "TMe-2196" "TMe-2791" "TMe-188"  "TMe-1053" "TMe-1174" "TMe-1769"
#> [37] "TMe-3095"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out5, use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "mcquitty")


## WPGMC
sel_hclust_medoid_out6 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid", hclust.method = "median")
sel_hclust_medoid_out6
#> $I
#>  [1] "TMe-1696" "TMe-41"   "TMe-44"   "TMe-117"  "TMe-410"  "TMe-469" 
#>  [7] "TMe-500"  "TMe-566"  "TMe-569"  "TMe-2604" "TMe-815"  "TMe-1096"
#> [13] "TMe-1190" "TMe-1425" "TMe-1451" "TMe-1466" "TMe-1564" "TMe-1922"
#> [19] "TMe-1935" "TMe-1958" "TMe-1973" "TMe-2010" "TMe-2027" "TMe-2372"
#> [25] "TMe-2773" "TMe-2967" "TMe-3151" "TMe-3272" "TMe-3345" "TMe-3353"
#> [31] "TMe-3452" "TMe-3481" "TMe-3633" "TMe-606"  "TMe-717"  "TMe-1218"
#> [37] "TMe-1589" "TMe-2912" "TMe-2927" "TMe-2934" "TMe-2937" "TMe-2985"
#> [43] "TMe-2996" "TMe-3026" "TMe-3112" "TMe-3163" "TMe-3292" "TMe-3389"
#> [49] "TMe-3415" "TMe-3433" "TMe-3437" "TMe-3457" "TMe-3471" "TMe-3475"
#> [55] "TMe-3480" "TMe-3514" "TMe-3705" "TMe-2910" "TMe-2943" "TMe-3462"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-1969" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-477"  "TMe-509" 
#>  [7] "TMe-768"  "TMe-1444" "TMe-1907" "TMe-1919" "TMe-2033" "TMe-2056"
#> [13] "TMe-2128" "TMe-2204" "TMe-2891" "TMe-2970" "TMe-3009" "TMe-3140"
#> [19] "TMe-3141" "TMe-3401" "TMe-3466" "TMe-3605" "TMe-47"   "TMe-1619"
#> [25] "TMe-2352" "TMe-2860" "TMe-2952" "TMe-2957" "TMe-3200" "TMe-3210"
#> [31] "TMe-3222" "TMe-3237" "TMe-3338" "TMe-3366" "TMe-3406" "TMe-3530"
#> [37] "TMe-3547" "TMe-3690" "TMe-3800" "TMe-3286" "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-2154" "TMe-104"  "TMe-141"  "TMe-261" 
#>  [7] "TMe-267"  "TMe-381"  "TMe-405"  "TMe-420"  "TMe-425"  "TMe-785" 
#> [13] "TMe-946"  "TMe-1200" "TMe-1868" "TMe-1897" "TMe-1993" "TMe-2802"
#> [19] "TMe-2823" "TMe-2897" "TMe-2968" "TMe-3234" "TMe-3299" "TMe-3804"
#> [25] "TMe-35"   "TMe-174"  "TMe-234"  "TMe-635"  "TMe-2756" "TMe-2926"
#> [31] "TMe-3143" "TMe-3274" "TMe-3383" "TMe-3397" "TMe-3445" "TMe-3565"
#> [37] "TMe-3575"
#> 
#> $IV
#>  [1] "TMe-2340" "TMe-154"  "TMe-266"  "TMe-274"  "TMe-427"  "TMe-540" 
#>  [7] "TMe-550"  "TMe-656"  "TMe-734"  "TMe-812"  "TMe-897"  "TMe-919" 
#> [13] "TMe-1162" "TMe-1350" "TMe-1376" "TMe-1406" "TMe-1434" "TMe-1525"
#> [19] "TMe-1776" "TMe-1806" "TMe-1820" "TMe-2025" "TMe-2043" "TMe-2755"
#> [25] "TMe-2809" "TMe-2956" "TMe-2971" "TMe-2979" "TMe-3025" "TMe-3054"
#> [31] "TMe-3055" "TMe-3089" "TMe-3144" "TMe-3196" "TMe-3255" "TMe-3256"
#> [37] "TMe-3273" "TMe-3428" "TMe-3545" "TMe-3576" "TMe-3752" "TMe-5"   
#> [43] "TMe-27"   "TMe-205"  "TMe-209"  "TMe-225"  "TMe-226"  "TMe-259" 
#> [49] "TMe-372"  "TMe-398"  "TMe-421"  "TMe-608"  "TMe-619"  "TMe-663" 
#> [55] "TMe-698"  "TMe-737"  "TMe-761"  "TMe-1010" "TMe-1020" "TMe-1078"
#> [61] "TMe-1129" "TMe-1511" "TMe-1597" "TMe-1988" "TMe-3097" "TMe-3189"
#> [67] "TMe-3297" "TMe-3386" "TMe-3422" "TMe-3527" "TMe-3558" "TMe-3619"
#> [73] "TMe-3701" "TMe-3707" "TMe-3727" "TMe-3730" "TMe-3053" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-332"  "TMe-294"  "TMe-1541" "TMe-362"  "TMe-373" 
#>  [7] "TMe-376"  "TMe-378"  "TMe-592"  "TMe-399"  "TMe-418"  "TMe-419" 
#> [13] "TMe-423"  "TMe-547"  "TMe-603"  "TMe-755"  "TMe-782"  "TMe-828" 
#> [19] "TMe-861"  "TMe-1099" "TMe-1204" "TMe-1366" "TMe-1375" "TMe-1381"
#> [25] "TMe-1427" "TMe-1458" "TMe-1482" "TMe-1963" "TMe-2041" "TMe-2050"
#> [31] "TMe-2060" "TMe-2121" "TMe-2531" "TMe-2589" "TMe-2688" "TMe-2790"
#> [37] "TMe-2853" "TMe-2855" "TMe-2863" "TMe-2973" "TMe-3185" "TMe-162" 
#> [43] "TMe-189"  "TMe-193"  "TMe-208"  "TMe-224"  "TMe-277"  "TMe-481" 
#> [49] "TMe-536"  "TMe-616"  "TMe-627"  "TMe-688"  "TMe-700"  "TMe-712" 
#> [55] "TMe-723"  "TMe-748"  "TMe-816"  "TMe-866"  "TMe-975"  "TMe-1042"
#> [61] "TMe-1196" "TMe-1283" "TMe-1332" "TMe-1453" "TMe-1785" "TMe-2003"
#> [67] "TMe-2139" "TMe-2355" "TMe-2905" "TMe-2953" "TMe-2959" "TMe-3034"
#> [73] "TMe-3311" "TMe-3356" "TMe-3408" "TMe-3628" "TMe-3736" "TMe-3659"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-662"  "TMe-269"  "TMe-315"  "TMe-442"  "TMe-531"  "TMe-1816"
#>  [7] "TMe-725"  "TMe-751"  "TMe-1383" "TMe-1403" "TMe-1416" "TMe-1428"
#> [13] "TMe-1445" "TMe-1518" "TMe-1633" "TMe-1775" "TMe-1875" "TMe-1995"
#> [19] "TMe-2026" "TMe-2035" "TMe-2067" "TMe-2534" "TMe-2791" "TMe-310" 
#> [25] "TMe-693"  "TMe-963"  "TMe-1053" "TMe-1174" "TMe-1232" "TMe-1558"
#> [31] "TMe-1560" "TMe-1769" "TMe-2510" "TMe-2543" "TMe-3095" "TMe-2983"
#> [37] "TMe-3116"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out6, use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "median")


## UPGMC
sel_hclust_medoid_out7 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid", hclust.method = "centroid")
sel_hclust_medoid_out7
#> $I
#>  [1] "TMe-1696" "TMe-41"   "TMe-132"  "TMe-299"  "TMe-410"  "TMe-500" 
#>  [7] "TMe-566"  "TMe-569"  "TMe-815"  "TMe-867"  "TMe-1360" "TMe-1425"
#> [13] "TMe-1466" "TMe-1564" "TMe-1786" "TMe-1922" "TMe-2027" "TMe-2069"
#> [19] "TMe-2372" "TMe-2773" "TMe-2785" "TMe-2945" "TMe-2967" "TMe-3030"
#> [25] "TMe-3151" "TMe-3169" "TMe-3272" "TMe-3330" "TMe-3452" "TMe-3481"
#> [31] "TMe-264"  "TMe-606"  "TMe-717"  "TMe-1218" "TMe-1532" "TMe-1589"
#> [37] "TMe-2927" "TMe-2934" "TMe-2984" "TMe-2985" "TMe-2993" "TMe-3112"
#> [43] "TMe-3163" "TMe-3282" "TMe-3292" "TMe-3389" "TMe-3415" "TMe-3433"
#> [49] "TMe-3437" "TMe-3457" "TMe-3471" "TMe-3475" "TMe-3485" "TMe-3488"
#> [55] "TMe-3705" "TMe-2940" "TMe-2943" "TMe-3249" "TMe-3462" "TMe-3478"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-3495" "TMe-160"  "TMe-289"  "TMe-369"  "TMe-509"  "TMe-539" 
#>  [7] "TMe-902"  "TMe-1385" "TMe-1444" "TMe-1907" "TMe-1919" "TMe-2033"
#> [13] "TMe-2054" "TMe-2056" "TMe-2128" "TMe-2204" "TMe-2891" "TMe-2970"
#> [19] "TMe-3009" "TMe-3140" "TMe-3141" "TMe-3466" "TMe-3605" "TMe-47"  
#> [25] "TMe-171"  "TMe-478"  "TMe-1619" "TMe-2352" "TMe-2568" "TMe-2952"
#> [31] "TMe-3200" "TMe-3210" "TMe-3222" "TMe-3258" "TMe-3338" "TMe-3366"
#> [37] "TMe-3530" "TMe-3690" "TMe-3771" "TMe-3800" "TMe-3805"
#> 
#> $III
#>  [1] "TMe-4"    "TMe-6"    "TMe-2154" "TMe-141"  "TMe-261"  "TMe-267" 
#>  [7] "TMe-381"  "TMe-405"  "TMe-425"  "TMe-785"  "TMe-946"  "TMe-1200"
#> [13] "TMe-1593" "TMe-1868" "TMe-1897" "TMe-1939" "TMe-2897" "TMe-2968"
#> [19] "TMe-3032" "TMe-3299" "TMe-3804" "TMe-35"   "TMe-116"  "TMe-174" 
#> [25] "TMe-635"  "TMe-878"  "TMe-2405" "TMe-2756" "TMe-2926" "TMe-3143"
#> [31] "TMe-3274" "TMe-3383" "TMe-3445" "TMe-3551" "TMe-3565" "TMe-3575"
#> [37] "TMe-3715"
#> 
#> $IV
#>  [1] "TMe-2340" "TMe-154"  "TMe-274"  "TMe-388"  "TMe-427"  "TMe-550" 
#>  [7] "TMe-650"  "TMe-656"  "TMe-766"  "TMe-812"  "TMe-897"  "TMe-919" 
#> [13] "TMe-962"  "TMe-1162" "TMe-1348" "TMe-1368" "TMe-1376" "TMe-1406"
#> [19] "TMe-1434" "TMe-1525" "TMe-1806" "TMe-2025" "TMe-2043" "TMe-2059"
#> [25] "TMe-2788" "TMe-2809" "TMe-2956" "TMe-2971" "TMe-2979" "TMe-3025"
#> [31] "TMe-3054" "TMe-3055" "TMe-3089" "TMe-3255" "TMe-3256" "TMe-3273"
#> [37] "TMe-3428" "TMe-3545" "TMe-3576" "TMe-3752" "TMe-27"   "TMe-180" 
#> [43] "TMe-186"  "TMe-205"  "TMe-209"  "TMe-226"  "TMe-259"  "TMe-280" 
#> [49] "TMe-372"  "TMe-398"  "TMe-421"  "TMe-608"  "TMe-619"  "TMe-663" 
#> [55] "TMe-698"  "TMe-737"  "TMe-1020" "TMe-1129" "TMe-1511" "TMe-1988"
#> [61] "TMe-2567" "TMe-2922" "TMe-3167" "TMe-3189" "TMe-3225" "TMe-3297"
#> [67] "TMe-3386" "TMe-3417" "TMe-3422" "TMe-3527" "TMe-3701" "TMe-3707"
#> [73] "TMe-3727" "TMe-3730" "TMe-2960" "TMe-3053" "TMe-3316" "TMe-3442"
#> [79] "TMe-3573" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-48"   "TMe-332"  "TMe-336"  "TMe-362"  "TMe-378"  "TMe-399" 
#>  [7] "TMe-418"  "TMe-419"  "TMe-423"  "TMe-547"  "TMe-603"  "TMe-782" 
#> [13] "TMe-861"  "TMe-972"  "TMe-1099" "TMe-1204" "TMe-1215" "TMe-1366"
#> [19] "TMe-1375" "TMe-1381" "TMe-1427" "TMe-1458" "TMe-1482" "TMe-1629"
#> [25] "TMe-1924" "TMe-1963" "TMe-2041" "TMe-2050" "TMe-2060" "TMe-2121"
#> [31] "TMe-2290" "TMe-2413" "TMe-2530" "TMe-2531" "TMe-2589" "TMe-2688"
#> [37] "TMe-2790" "TMe-2853" "TMe-2863" "TMe-2973" "TMe-3185" "TMe-189" 
#> [43] "TMe-208"  "TMe-284"  "TMe-481"  "TMe-536"  "TMe-616"  "TMe-655" 
#> [49] "TMe-700"  "TMe-707"  "TMe-712"  "TMe-723"  "TMe-748"  "TMe-759" 
#> [55] "TMe-798"  "TMe-803"  "TMe-816"  "TMe-866"  "TMe-877"  "TMe-975" 
#> [61] "TMe-1042" "TMe-1196" "TMe-1199" "TMe-1305" "TMe-1311" "TMe-1332"
#> [67] "TMe-2003" "TMe-2139" "TMe-2355" "TMe-2905" "TMe-2953" "TMe-2959"
#> [73] "TMe-3034" "TMe-3311" "TMe-3356" "TMe-3408" "TMe-3736" "TMe-3659"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-130"  "TMe-1550" "TMe-269"  "TMe-361"  "TMe-442"  "TMe-465" 
#>  [7] "TMe-725"  "TMe-751"  "TMe-1383" "TMe-1403" "TMe-1416" "TMe-1428"
#> [13] "TMe-1445" "TMe-1518" "TMe-1875" "TMe-2035" "TMe-2064" "TMe-2067"
#> [19] "TMe-2196" "TMe-2534" "TMe-2791" "TMe-236"  "TMe-310"  "TMe-693" 
#> [25] "TMe-1053" "TMe-1124" "TMe-1174" "TMe-1182" "TMe-1232" "TMe-1769"
#> [31] "TMe-2510" "TMe-2543" "TMe-3095" "TMe-3387" "TMe-2983" "TMe-3116"
#> [37] "TMe-2551"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out7, use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "centroid")


## Ward's D2
sel_hclust_medoid_out8 <-
  select.distance(data = data, names = "genotypes",
                  group = "Cluster", alloc = counts, dist.mat = dist_matrix,
                  always.selected = mand_accns,
                  method = "hclust.medoid", hclust.method = "ward.D2")
sel_hclust_medoid_out8
#> $I
#>  [1] "TMe-642"  "TMe-41"   "TMe-3149" "TMe-300"  "TMe-3424" "TMe-2131"
#>  [7] "TMe-1564" "TMe-1931" "TMe-3202" "TMe-489"  "TMe-500"  "TMe-2255"
#> [13] "TMe-3441" "TMe-1448" "TMe-1981" "TMe-1341" "TMe-756"  "TMe-2253"
#> [19] "TMe-1734" "TMe-3132" "TMe-3330" "TMe-1823" "TMe-610"  "TMe-882" 
#> [25] "TMe-1170" "TMe-910"  "TMe-1086" "TMe-1964" "TMe-1955" "TMe-306" 
#> [31] "TMe-1621" "TMe-1717" "TMe-1940" "TMe-3223" "TMe-2810" "TMe-3477"
#> [37] "TMe-3490" "TMe-3729" "TMe-3548" "TMe-3296" "TMe-3433" "TMe-3220"
#> [43] "TMe-3268" "TMe-2985" "TMe-2916" "TMe-2917" "TMe-3479" "TMe-3236"
#> [49] "TMe-3685" "TMe-248"  "TMe-2943" "TMe-3501" "TMe-2927" "TMe-3493"
#> [55] "TMe-2996" "TMe-3419" "TMe-3398" "TMe-3252" "TMe-3425" "TMe-3460"
#> [61] "TMe-3423"
#> 
#> $II
#>  [1] "TMe-40"   "TMe-2611" "TMe-3571" "TMe-85"   "TMe-3009" "TMe-2414"
#>  [7] "TMe-377"  "TMe-2715" "TMe-648"  "TMe-2765" "TMe-1474" "TMe-455" 
#> [13] "TMe-2329" "TMe-706"  "TMe-2951" "TMe-589"  "TMe-601"  "TMe-2843"
#> [19] "TMe-2824" "TMe-2200" "TMe-674"  "TMe-2123" "TMe-2258" "TMe-1698"
#> [25] "TMe-176"  "TMe-1833" "TMe-1505" "TMe-2021" "TMe-2568" "TMe-2033"
#> [31] "TMe-3498" "TMe-2995" "TMe-74"   "TMe-3455" "TMe-3642" "TMe-3093"
#> [37] "TMe-3540" "TMe-3766" "TMe-3605" "TMe-2952" "TMe-3188"
#> 
#> $III
#>  [1] "TMe-161"  "TMe-3029" "TMe-13"   "TMe-1836" "TMe-64"   "TMe-3048"
#>  [7] "TMe-187"  "TMe-200"  "TMe-206"  "TMe-1138" "TMe-3005" "TMe-267" 
#> [13] "TMe-420"  "TMe-3028" "TMe-853"  "TMe-2394" "TMe-946"  "TMe-1790"
#> [19] "TMe-1398" "TMe-3118" "TMe-3277" "TMe-3592" "TMe-1207" "TMe-3721"
#> [25] "TMe-2158" "TMe-2356" "TMe-1288" "TMe-3094" "TMe-3324" "TMe-3207"
#> [31] "TMe-3750" "TMe-3043" "TMe-3128" "TMe-3458" "TMe-23"   "TMe-3317"
#> [37] "TMe-3143"
#> 
#> $IV
#>  [1] "TMe-3073" "TMe-3045" "TMe-204"  "TMe-3084" "TMe-3357" "TMe-72"  
#>  [7] "TMe-525"  "TMe-1364" "TMe-3283" "TMe-207"  "TMe-318"  "TMe-270" 
#> [13] "TMe-594"  "TMe-314"  "TMe-3204" "TMe-82"   "TMe-1406" "TMe-386" 
#> [19] "TMe-387"  "TMe-1148" "TMe-396"  "TMe-1726" "TMe-428"  "TMe-3168"
#> [25] "TMe-1579" "TMe-513"  "TMe-3435" "TMe-87"   "TMe-550"  "TMe-3542"
#> [31] "TMe-840"  "TMe-1700" "TMe-1430" "TMe-3757" "TMe-3019" "TMe-834" 
#> [37] "TMe-1397" "TMe-18"   "TMe-1139" "TMe-1402" "TMe-2039" "TMe-3727"
#> [43] "TMe-1151" "TMe-1419" "TMe-802"  "TMe-1342" "TMe-2399" "TMe-1350"
#> [49] "TMe-190"  "TMe-3409" "TMe-885"  "TMe-1027" "TMe-1485" "TMe-1336"
#> [55] "TMe-3537" "TMe-2332" "TMe-1903" "TMe-2928" "TMe-2958" "TMe-2809"
#> [61] "TMe-210"  "TMe-3191" "TMe-2979" "TMe-3231" "TMe-3615" "TMe-3068"
#> [67] "TMe-3594" "TMe-3535" "TMe-3108" "TMe-3562" "TMe-3255" "TMe-3253"
#> [73] "TMe-3581" "TMe-3044" "TMe-2947" "TMe-153"  "TMe-2758" "TMe-699" 
#> [79] "TMe-1278" "TMe-34"   "TMe-801" 
#> 
#> $V
#>  [1] "TMe-3332" "TMe-99"   "TMe-1127" "TMe-167"  "TMe-212"  "TMe-1339"
#>  [7] "TMe-294"  "TMe-305"  "TMe-945"  "TMe-1399" "TMe-1979" "TMe-1358"
#> [13] "TMe-2441" "TMe-826"  "TMe-474"  "TMe-982"  "TMe-1067" "TMe-1099"
#> [19] "TMe-629"  "TMe-385"  "TMe-390"  "TMe-1011" "TMe-417"  "TMe-685" 
#> [25] "TMe-217"  "TMe-1458" "TMe-1859" "TMe-430"  "TMe-1609" "TMe-476" 
#> [31] "TMe-487"  "TMe-870"  "TMe-2845" "TMe-585"  "TMe-884"  "TMe-623" 
#> [37] "TMe-645"  "TMe-1312" "TMe-730"  "TMe-341"  "TMe-2355" "TMe-844" 
#> [43] "TMe-1877" "TMe-2213" "TMe-2835" "TMe-813"  "TMe-651"  "TMe-673" 
#> [49] "TMe-1523" "TMe-1215" "TMe-1188" "TMe-2933" "TMe-1307" "TMe-1901"
#> [55] "TMe-1234" "TMe-997"  "TMe-1440" "TMe-1482" "TMe-1500" "TMe-2753"
#> [61] "TMe-1788" "TMe-2192" "TMe-858"  "TMe-1269" "TMe-2050" "TMe-2326"
#> [67] "TMe-850"  "TMe-2530" "TMe-2290" "TMe-1487" "TMe-2435" "TMe-3185"
#> [73] "TMe-669"  "TMe-284"  "TMe-1227" "TMe-771"  "TMe-1733" "TMe-2905"
#> [79] "TMe-2018" "TMe-551" 
#> 
#> $VI
#>  [1] "TMe-683"  "TMe-1985" "TMe-678"  "TMe-1477" "TMe-1945" "TMe-620" 
#>  [7] "TMe-2308" "TMe-626"  "TMe-548"  "TMe-1110" "TMe-2402" "TMe-273" 
#> [13] "TMe-281"  "TMe-1481" "TMe-1164" "TMe-846"  "TMe-876"  "TMe-1079"
#> [19] "TMe-1460" "TMe-2963" "TMe-625"  "TMe-1413" "TMe-1416" "TMe-1661"
#> [25] "TMe-1445" "TMe-1486" "TMe-1518" "TMe-2389" "TMe-1775" "TMe-1832"
#> [31] "TMe-2551" "TMe-2035" "TMe-1608" "TMe-2196" "TMe-188"  "TMe-1053"
#> [37] "TMe-3095"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_hclust_medoid_out8, use.names = FALSE)) +
  labs(title = "hclust.medoid", subtitle = "ward.D2")
```
