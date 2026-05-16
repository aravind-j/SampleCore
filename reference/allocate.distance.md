# Allocation of Entries to be Selected from Clusters/Groups based on Distance-based Diversity Metrics for Core Collection Development

Estimate the number of entries to be allocated from each cluster/group
in the entire collection to construct a core collection on the basis of
different metrics computed from within cluster/group distances. The
following strategies are implemented.

- Diversity (Distance based)

- Diversity (Distance based) & Proportional

- Diversity (Distance based) & Logarithmic

- Diversity (Distance based) & Square root

## Usage

``` r
allocate.distance(
  data,
  names,
  group,
  dist.mat,
  method = c("dist", "dist.prop", "dist.log", "dist.sqrt"),
  metric = c("mean", "median", "max", "range", "mnnd", "mdc", "mdm", "mstl", "nclust"),
  clust.fun = NULL,
  log.base = exp(1),
  size
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

- dist.mat:

  A precomputed distance matrix of distance measures between the
  accessions in `data`.

- method:

  The allocation method. Either `"dist"` for constant or `"dist.prop"`
  for proportional or `"dist.log"` for logarithmic or `"dist.sqrt"` for
  square root allocation. See **Methods**.

- metric:

  The metric to be computed from the distance matrix. Either `"mean"`,
  `"median"`, `"max"`, `"range"`, `"mnnd"`, `"mdc"`, `"mdm"`, `"mstl"`,
  or `"nclust"`. See **Metrics**.

- clust.fun:

  A function to generate clusters from a distance matrix and return the
  number of clusters.

- log.base:

  The logarithm base to be used for logarithmic method of sampling.
  Default is `exp(1)`.

- size:

  The desired core set size proportion.

## Value

A named numeric vector specifying the number of entries to be selected
from each cluster/group. The vector names correspond to the levels of
the "`"group"` column, and values indicate the number of elements to be
selected from each level.

## Details

The number of entries to be chosen from each cluster is estimated either
on the basis of diversity of entries within that cluster/group alone or
in combination with the size of the cluster/group (See **Methods**).

The within-cluster/group diversity is estimated as several metrics from
the within cluster/group genetic distances between accessions (See
**Metrics**).

Franco et al. (2005) proposed a method based on mean Gower's distance
(Gower 1971) which was also extended to other distance measure averages
named D Allocation strategy (Franco et al. 2006) . These methods were
also combined with the proportional and logarithmic methods. For
example, the GP and GL strategy of Bisht et al. (1999) and Mahajan et
al. (1999) as well as the NY and LD allocation methods of Franco et al.
(2005) .

## Methods

### Diversity method

From an entire collection of size \\N\\, to construct a core set of
sample size \\n\\, the number of entries to be selected from the \\i\\th
group among \\1 \cdots g\\ groups (\\n\_{i}\\) is estimated as below.

\\n\_{i} = n \times \frac{D\_{i}}{\sum\_{i=1}^{g}D\_{i}}\\

Where, \\D\_{i}\\ is a measure of the extent of diversity present in the
\\i\\th cluster.

### Diversity and proportional method

Here the number of entries to be selected is proportional to the
diversity of the cluster/group (\\D\_{i}\\) weighted by the the
cluster/group size (\\N\_{i}\\).

\\n\_{i} = n \times \frac{N\_{i}D\_{i}}{\sum\_{i=1}^{g}N\_{i}D\_{i}}\\

### Diversity and logarithmic method

Here the number of entries to be selected is proportional to the
diversity of the cluster/group (\\D\_{i}\\) weighted by the logarithm of
the cluster/group size (\\N\_{i}\\).

\\n\_{i} = n \times
\frac{\log(N\_{i})D\_{i}}{\sum\_{i=1}^{g}\log(N\_{i})D\_{i}}\\

### Diversity and square root method

Here the number of entries to be selected is proportional to the
diversity of the cluster/group (\\D\_{i}\\) weighted by the square root
of the cluster/group size (\\N\_{i}\\).

\\n\_{i} = n \times
\frac{\sqrt{N\_{i}}D\_{i}}{\sum\_{i=1}^{g}\sqrt{N\_{i}}D\_{i}}\\

## Metrics

### Summary/Decriptive statistics

These include mean, median, maximum and range of genetic distances
between entries in a cluster.

### Mean nearest-neighbour distance (\\MNND\\)

It is the average, across all entries, of the distance to each entry’s
closest other entry (\\d\_{g\_{min}}\\), based on a genetic given
distance matrix (Clark and Evans 1954) .

For each entry, the nearest-neighbour distance (\\d\_{g\_{min}}\\) is
the smallest non-zero distance with any other entry.

\\d\_{g\_{min}} = \min\_{h \ne g} d\_{gh}\\

The Mean nearest-neighbour distance (\\MNND\\) can then be computed as:

\\\textrm{MNND} = \frac{1}{G} \sum\_{g=1}^{G} d_g\\

Where, (\\g\\) is the index of an entry in a genetic distance matrix,
\\h\\ is the index of all other genotypes and \\G\\ is the total number
of genotypes in a cluster/group.

### Minimum spanning tree length (\\MSTL\\)

It is defined as the sum of edge weights in the minimum spanning tree
constructed from the genetic distance matrix of entries within a
cluster/group. A minimum spanning tree (MST) connects all entries such
that the total distance is minimized and no cycles are formed. It
represents the most efficient way to connect all entries based on
pairwise genetic distances (Gower and Ross 1969) .

For genetic distance \\d\_{gh}\\ between entries \\g\\ and \\h\\, the
MST is a subset of edges that connects all \\G\\ entries with exactly
\\G - 1\\ edges and minimum total weight. The MST length (\\MSTL\\) can
then be computed as:

\\\textrm{MSTL} = \sum\_{(g,h) \in \mathcal{T}} d\_{gh}\\

Where \\\mathcal{T}\\ denotes the set of edges in the MST.

### Mean distance to centroid and median (\\MDC\\, \\MDM\\)

These quantify the average dispersion of entries within a cluster/group
relative to a central point in multivariate space derived from the
genetic distance matrix.

The centroid represents the multivariate mean position of all entries in
a cluster (Sokal and Sneath 1963; Sneath and Sokal 1973) ., whereas the
median (spatial median) provides a robust central location that is less
influenced by extreme values (Bradley et al. 1999) .

For \\d\_{gC}\\ and \\d\_{gM}\\ distances of entry \\g\\ from the
centroid \\C\\ and median \\M\\, respectively. These measures are
computed as:

\\\textrm{MDC} = \frac{1}{G} \sum\_{g=1}^{G} d\_{gC}\\

\\\textrm{MDM} = \frac{1}{G} \sum\_{g=1}^{G} d\_{gM}\\

Where \\G\\ is the total number of entries in the cluster/group.

### Number of clusters

(Diwan et al. 1994) proposed the number of clusters produced by a
multivariate cluster analysis at a specific distance threshold as an
estimate of the diversity.

## References

Bisht IS, Mahajan RK, Gautam PL (1999). “Assessment of genetic
diversity, stratification of germplasm accessions in diversity groups
and sampling strategies for establishing a core collection of Indian
sesame (*Sesamum indicum* L.).” *Plant Genetic Resources Newsletter*,
**199 Supp.**, 35–46.  
  
Bradley PS, Bennett KP, Mangasarian OL (1999). “Constrained k-means
clustering.” Technical Report MSR-TR-2000-65, Microsoft Research,
Redmond, WA.  
  
Clark PJ, Evans FC (1954). “Distance to nearest neighbor as a measure of
spatial relationships in populations.” *Ecology*, **35**(4), 445–453.  
  
Diwan N, Bauchan GR, McIntosh MS (1994). “A core collection for the
united states annual *Medicago* germplasm collection.” *Crop Science*,
**34**(1), cropsci1994.0011183X003400010051x.  
  
Franco J, Crossa J, Taba S, Shands H (2005). “A sampling strategy for
conserving genetic diversity when forming core subsets.” *Crop Science*,
**45**(3), 1035–1044.  
  
Franco J, Crossa J, Warburton ML, Taba S (2006). “Sampling strategies
for conserving maize diversity when forming core subsets using genetic
markers.” *Crop Science*, **46**(2), 854–864.  
  
Gower JC (1971). “A general coefficient of similarity and some of its
properties.” *Biometrics*, **27**(4), 857–871.  
  
Gower JC, Ross GJS (1969). “Minimum spanning trees and single linkage
cluster analysis.” *Journal of the Royal Statistical Society. Series C
(Applied Statistics)*, **18**(1), 54–64.  
  
Mahajan RK, Bisht IS, Gautam PL (1999). “Sampling strategies for
developing Indian sesame core collection.” *Indian Journal of Plant
Genetic Resources*, **12**(01), 1–9.  
  
Sneath PHA, Sokal RR (1973). *Numerical Taxonomy: The Principles and
Practice of Numerical Classification*, A Series of books in biology. W.
H. Freeman, San Francisco. ISBN 978-0-7167-0697-7.  
  
Sokal RR, Sneath PHA (1963). *Principles of numerical taxonomy*, A
Series of books in biology. W. H. Freeman, San Francisco.

## See also

[`allocate.basic`](https://aravind-j.github.io/SampleCore/reference/allocate.basic.md),
[`allocate.diversity`](https://aravind-j.github.io/SampleCore/reference/allocate.diversity.md)

## Examples

``` r
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare example data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(cluster)

# Get distance matrix
data("cassava_EC_gp")

set.seed(123)
cassava_EC_gp <- cassava_EC_gp[sample(1:nrow(cassava_EC_gp), 500), ]

quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW",
           "AVPW", "ARSR", "SRDM")
qual <- c("CUAL", "LNGS", "PTLC", "DSTA", "LFRT", "LBTEF", "CBTR", "NMLB",
          "ANGB", "CUAL9M", "LVC9M", "TNPR9M", "PL9M", "STRP", "STRC",
          "PSTR")

data <- cassava_EC_gp

# Convert qualitative data columns to factor
data[, qual] <- lapply(data[, qual], as.factor)

# Standardise quantitative data column
data[, quant] <- lapply(data[, quant], function(x) {
  scale(x)[, 1]
})

# Get the Gower's distance matrix
dist_matrix <- daisy(x = data[, c(qual, quant)],
                     metric = "gower")

# Get data
data <- cassava_EC_gp
data <- cbind(genotypes = rownames(cassava_EC_gp), cassava_EC_gp)
row.names(data) <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Custom clustering functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# UPGMA with hclust
clust_fun_upgma <- function(x) {
  # Tree
  tree_out <- hclust(x, method = "average")
  # Clusters
  cutree(tree_out, h = 0.2)
}

if (requireNamespace('fastcluster', quietly = TRUE)) {
  # Ward's minimum variance with fastcluster
  clust_fun_ward <- function(x) {
    # Tree
    tree_out <- fastcluster::hclust(x, method = "ward.D2")
    # Clusters
    cutree(tree_out, h = 0.2)
  }
}

if (requireNamespace('dbscan', quietly = TRUE)) {
  # Density-based clustering with dbscan
  clust_fun_dbscan <- function(x) {
    clust_out <- dbscan::dbscan(x, eps = 0.25)
    # remove noise: TODO
    setNames(clust_out$cluster, labels(x))
  }
}

if (requireNamespace('biotools', quietly = TRUE)) {
  # Tocher's sequential clustering
  clust_fun_tocher <- function(x) {
    clust_out <- biotools::tocher(x, algorithm = "sequential")
    setNames(clust_out$class, labels(x))
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity (Distance based) allocation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Mean
dist_out_mean <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist",
                    dist.mat = dist_matrix, metric = "mean",
                    size = 0.2)
dist_out_mean
#>   I  II III  IV   V  VI 
#>  18  13  16  18  15  20 

## Median
dist_out_median <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist",
                    dist.mat = dist_matrix, metric = "median",
                    size = 0.2)
dist_out_median
#>   I  II III  IV   V  VI 
#>  18  13  16  18  15  20 

## Maximum
dist_out_max <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist",
                    dist.mat = dist_matrix, metric = "max",
                    size = 0.2)
dist_out_max
#>   I  II III  IV   V  VI 
#>  18  12  18  17  15  20 

## Range
dist_out_range <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist",
                    dist.mat = dist_matrix, metric = "range",
                    size = 0.2)
dist_out_range
#>   I  II III  IV   V  VI 
#>  17  12  18  19  15  20 

## Mean nearest-neighbour distance
dist_out_mnnd <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist",
                    dist.mat = dist_matrix, metric = "mnnd",
                    size = 0.2)
dist_out_mnnd
#>   I  II III  IV   V  VI 
#>  20  14  15  17  15  20 

## Minimum spanning tree length
dist_out_mstl <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist",
                    dist.mat = dist_matrix, metric = "mstl",
                    size = 0.2)
dist_out_mstl
#>   I  II III  IV   V  VI 
#>  11   9  21  28  10  21 

# \donttest{
  ## Mean distance to centroid
  dist_out_mdc <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist",
                      dist.mat = dist_matrix, metric = "mdc",
                      size = 0.2)
  dist_out_mdc
#>   I  II III  IV   V  VI 
#>  20  15  13  18  16  18 

  ## Mean distance to median
  dist_out_mdm <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist",
                      dist.mat = dist_matrix, metric = "mdm",
                      size = 0.2)
  dist_out_mdm
#>   I  II III  IV   V  VI 
#>  20  15  13  18  16  18 
# }

## Number of clusters

### UPGMA with hclust
dist_out_nclust1 <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist",
                    dist.mat = dist_matrix, metric = "nclust",
                    clust.fun = clust_fun_upgma,
                    size = 0.2)
dist_out_nclust1
#>   I  II III  IV   V  VI 
#>  12   7  20  28  11  22 

# Ward's minimum variance with fastcluster
if (requireNamespace('fastcluster', quietly = TRUE)) {
  dist_out_nclust2 <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist",
                      dist.mat = dist_matrix, metric = "nclust",
                      clust.fun = clust_fun_ward,
                      size = 0.2)
  dist_out_nclust2
}
#>   I  II III  IV   V  VI 
#>  12   8  21  27  11  21 


# Density-based clustering with dbscan
if (requireNamespace('dbscan', quietly = TRUE)) {
  dist_out_nclust3 <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist",
                      dist.mat = dist_matrix, metric = "nclust",
                      clust.fun = clust_fun_dbscan,
                      size = 0.2)
  dist_out_nclust3
}
#>   I  II III  IV   V  VI 
#>  18   9  18  18   9  27 

# \donttest{
  if (requireNamespace('biotools', quietly = TRUE)) {
    # Tocher's sequential clustering
    dist_out_nclust4 <-
      allocate.distance(data = data, names = "genotypes",
                        group = "Cluster", method = "dist",
                        dist.mat = dist_matrix, metric = "nclust",
                        clust.fun = clust_fun_tocher,
                        size = 0.2)
    dist_out_nclust4
  }
#>   I  II III  IV   V  VI 
#>  13  20  17  22  13  15 
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity (Distance based) & Proportional
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Mean
dist_prop_out_mean <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.prop",
                    dist.mat = dist_matrix, metric = "mean",
                    size = 0.2)
dist_prop_out_mean
#>   I  II III  IV   V  VI 
#>  18   9  11  29  22  11 

## Median
dist_prop_out_median <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.prop",
                    dist.mat = dist_matrix, metric = "median",
                    size = 0.2)
dist_prop_out_median
#>   I  II III  IV   V  VI 
#>  18   9  11  29  22  11 

## Maximum
dist_prop_out_max <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.prop",
                    dist.mat = dist_matrix, metric = "max",
                    size = 0.2)
dist_prop_out_max
#>   I  II III  IV   V  VI 
#>  19   8  12  29  21  11 

## Range
dist_prop_out_range <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.prop",
                    dist.mat = dist_matrix, metric = "range",
                    size = 0.2)
dist_prop_out_range
#>   I  II III  IV   V  VI 
#>  17   8  12  31  21  11 

## Mean nearest-neighbour distance
dist_prop_out_mnnd <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.prop",
                    dist.mat = dist_matrix, metric = "mnnd",
                    size = 0.2)
dist_prop_out_mnnd
#>   I  II III  IV   V  VI 
#>  20   9  10  28  21  12 

## Minimum spanning tree length
dist_prop_out_mstl <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.prop",
                    dist.mat = dist_matrix, metric = "mstl",
                    size = 0.2)
dist_prop_out_mstl
#>   I  II III  IV   V  VI 
#>  11   6  14  44  14  12 

# \donttest{
  ## Mean distance to centroid
  dist_prop_out_mdc <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.prop",
                      dist.mat = dist_matrix, metric = "mdc",
                      size = 0.2)
  dist_prop_out_mdc
#>   I  II III  IV   V  VI 
#>  20  10   9  29  22  10 

  ## Mean distance to median
  dist_prop_out_mdm <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.prop",
                      dist.mat = dist_matrix, metric = "mdm",
                      size = 0.2)
  dist_prop_out_mdm
#>   I  II III  IV   V  VI 
#>  20  10   9  29  22  10 
# }

## Number of clusters

### UPGMA with hclust
dist_prop_out_nclust1 <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.prop",
                    dist.mat = dist_matrix, metric = "nclust",
                    clust.fun = clust_fun_upgma,
                    size = 0.2)
dist_prop_out_nclust1
#>   I  II III  IV   V  VI 
#>  12   5  13  44  14  12 

# Ward's minimum variance with fastcluster
if (requireNamespace('fastcluster', quietly = TRUE)) {
  dist_prop_out_nclust2 <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.prop",
                      dist.mat = dist_matrix, metric = "nclust",
                      clust.fun = clust_fun_ward,
                      size = 0.2)
  dist_prop_out_nclust2
}
#>   I  II III  IV   V  VI 
#>  12   6  14  43  15  11 

# Density-based clustering with dbscan
if (requireNamespace('dbscan', quietly = TRUE)) {
  dist_prop_out_nclust3 <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.prop",
                      dist.mat = dist_matrix, metric = "nclust",
                      clust.fun = clust_fun_dbscan,
                      size = 0.2)
  dist_prop_out_nclust3
}
#>   I  II III  IV   V  VI 
#>  20   7  13  31  13  16 

# \donttest{
  if (requireNamespace('biotools', quietly = TRUE)) {
    # Tocher's sequential clustering
    dist_prop_out_nclust4 <-
      allocate.distance(data = data, names = "genotypes",
                        group = "Cluster", method = "dist.prop",
                        dist.mat = dist_matrix, metric = "nclust",
                        clust.fun = clust_fun_tocher,
                        size = 0.2)
    dist_prop_out_nclust4
  }
#>   I  II III  IV   V  VI 
#>  14  13  11  35  18   8 
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity (Distance based) & Logarithmic
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Mean
dist_log_out_mean <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.log",
                    dist.mat = dist_matrix, metric = "mean",
                    size = 0.2)
dist_log_out_mean
#>   I  II III  IV   V  VI 
#>  18  12  15  20  17  18 

## Median
dist_log_out_median <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.log",
                    dist.mat = dist_matrix, metric = "median",
                    size = 0.2)
dist_log_out_median
#>   I  II III  IV   V  VI 
#>  18  12  15  20  17  18 

## Maximum
dist_log_out_max <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.log",
                    dist.mat = dist_matrix, metric = "max",
                    size = 0.2)
dist_log_out_max
#>   I  II III  IV   V  VI 
#>  18  11  16  20  16  18 

## Range
dist_log_out_range <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.log",
                    dist.mat = dist_matrix, metric = "range",
                    size = 0.2)
dist_log_out_range
#>   I  II III  IV   V  VI 
#>  17  11  16  21  17  18 

## Mean nearest-neighbour distance
dist_log_out_mnnd <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.log",
                    dist.mat = dist_matrix, metric = "mnnd",
                    size = 0.2)
dist_log_out_mnnd
#>   I  II III  IV   V  VI 
#>  20  13  14  19  16  18 

## Minimum spanning tree length
dist_log_out_mstl <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.log",
                    dist.mat = dist_matrix, metric = "mstl",
                    size = 0.2)
dist_log_out_mstl
#>   I  II III  IV   V  VI 
#>  11   8  19  31  11  19 

# \donttest{
  ## Mean distance to centroid
  dist_log_out_mdc <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.log",
                      dist.mat = dist_matrix, metric = "mdc",
                      size = 0.2)
  dist_log_out_mdc
#>   I  II III  IV   V  VI 
#>  20  14  12  20  18  15 

  ## Mean distance to median
  dist_log_out_mdm <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.log",
                      dist.mat = dist_matrix, metric = "mdm",
                      size = 0.2)
  dist_log_out_mdm
#>   I  II III  IV   V  VI 
#>  20  14  12  20  18  15 
# }

## Number of clusters

### UPGMA with hclust
dist_log_out_nclust1 <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.log",
                    dist.mat = dist_matrix, metric = "nclust",
                    clust.fun = clust_fun_upgma,
                    size = 0.2)
dist_log_out_nclust1
#>   I  II III  IV   V  VI 
#>  13   7  18  31  12  20 

# Ward's minimum variance with fastcluster
if (requireNamespace('fastcluster', quietly = TRUE)) {
  dist_log_out_nclust2 <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.log",
                      dist.mat = dist_matrix, metric = "nclust",
                      clust.fun = clust_fun_ward,
                      size = 0.2)
  dist_log_out_nclust2
}
#>   I  II III  IV   V  VI 
#>  12   8  20  30  12  18 

if (requireNamespace('dbscan', quietly = TRUE)) {
  # Density-based clustering with dbscan
  dist_log_out_nclust3 <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.log",
                      dist.mat = dist_matrix, metric = "nclust",
                      clust.fun = clust_fun_dbscan,
                      size = 0.2)
  dist_log_out_nclust3
}
#>   I  II III  IV   V  VI 
#>  19   9  17  21  10  24 

# \donttest{
  if (requireNamespace('biotools', quietly = TRUE)) {
    # Tocher's sequential clustering
    dist_log_out_nclust4 <-
      allocate.distance(data = data, names = "genotypes",
                        group = "Cluster", method = "dist.log",
                        dist.mat = dist_matrix, metric = "nclust",
                        clust.fun = clust_fun_tocher,
                        size = 0.2)
    dist_log_out_nclust4
  }
#>   I  II III  IV   V  VI 
#>  14  19  15  24  15  13 
# }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity (Distance based) & Square root
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Mean
dist_sqrt_out_mean <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.sqrt",
                    dist.mat = dist_matrix, metric = "mean",
                    size = 0.2)
dist_sqrt_out_mean
#>   I  II III  IV   V  VI 
#>  18  11  13  23  19  15 

## Median
dist_sqrt_out_median <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.sqrt",
                    dist.mat = dist_matrix, metric = "median",
                    size = 0.2)
dist_sqrt_out_median
#>   I  II III  IV   V  VI 
#>  18  11  14  23  19  15 

## Maximum
dist_sqrt_out_max <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.sqrt",
                    dist.mat = dist_matrix, metric = "max",
                    size = 0.2)
dist_sqrt_out_max
#>   I  II III  IV   V  VI 
#>  19  10  15  23  18  15 

## Range
dist_sqrt_out_range <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.sqrt",
                    dist.mat = dist_matrix, metric = "range",
                    size = 0.2)
dist_sqrt_out_range
#>   I  II III  IV   V  VI 
#>  17  10  15  24  19  15 

## Mean nearest-neighbour distance
dist_sqrt_out_mnnd <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.sqrt",
                    dist.mat = dist_matrix, metric = "mnnd",
                    size = 0.2)
dist_sqrt_out_mnnd
#>   I  II III  IV   V  VI 
#>  20  11  13  22  18  16 

## Minimum spanning tree length
dist_sqrt_out_mstl <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.sqrt",
                    dist.mat = dist_matrix, metric = "mstl",
                    size = 0.2)
dist_sqrt_out_mstl
#>   I  II III  IV   V  VI 
#>  11   8  17  36  12  16 

# \donttest{
  ## Mean distance to centroid
  dist_sqrt_out_mdc <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.sqrt",
                      dist.mat = dist_matrix, metric = "mdc",
                      size = 0.2)
  dist_sqrt_out_mdc
#>   I  II III  IV   V  VI 
#>  20  13  11  23  19  13 

  ## Mean distance to median
  dist_sqrt_out_mdm <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.sqrt",
                      dist.mat = dist_matrix, metric = "mdm",
                      size = 0.2)
  dist_sqrt_out_mdm
#>   I  II III  IV   V  VI 
#>  20  13  11  23  19  13 
# }

## Number of clusters

### UPGMA with hclust
dist_sqrt_out_nclust1 <-
  allocate.distance(data = data, names = "genotypes",
                    group = "Cluster", method = "dist.sqrt",
                    dist.mat = dist_matrix, metric = "nclust",
                    clust.fun = clust_fun_upgma,
                    size = 0.2)
dist_sqrt_out_nclust1
#>   I  II III  IV   V  VI 
#>  13   6  16  36  13  17 

# Ward's minimum variance with fastcluster
if (requireNamespace('fastcluster', quietly = TRUE)) {
  dist_sqrt_out_nclust2 <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.sqrt",
                      dist.mat = dist_matrix, metric = "nclust",
                      clust.fun = clust_fun_ward,
                      size = 0.2)
  dist_sqrt_out_nclust2
}
#>   I  II III  IV   V  VI 
#>  12   7  18  35  13  16 

if (requireNamespace('dbscan', quietly = TRUE)) {
  # Density-based clustering with dbscan
  dist_sqrt_out_nclust3 <-
    allocate.distance(data = data, names = "genotypes",
                      group = "Cluster", method = "dist.sqrt",
                      dist.mat = dist_matrix, metric = "nclust",
                      clust.fun = clust_fun_dbscan,
                      size = 0.2)
  dist_sqrt_out_nclust3
}
#>   I  II III  IV   V  VI 
#>  19   8  16  24  11  21 

# \donttest{
  if (requireNamespace('biotools', quietly = TRUE)) {
    # Tocher's sequential clustering
    dist_sqrt_out_nclust4 <-
      allocate.distance(data = data, names = "genotypes",
                        group = "Cluster", method = "dist.sqrt",
                        dist.mat = dist_matrix, metric = "nclust",
                        clust.fun = clust_fun_tocher,
                        size = 0.2)
    dist_sqrt_out_nclust4
  }
#>   I  II III  IV   V  VI 
#>  14  17  14  28  16  11 
# }
```
