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
