# Allocation of Entries to be Selected from Clusters/Groups based on Diversity Index Estimates for Core Collection Development

Estimate the number of entries to be allocated from each cluster/group
in the entire collection to construct a core collection on the basis of
different metrics computed from within cluster/group diversity index
estimates. The following strategies are implemented.

- Diversity

- Diversity & Proportional

- Diversity & Logarithmic

- Diversity & Square root

## Usage

``` r
allocate.diversity(
  data,
  names,
  group,
  qualitative,
  method = c("div", "div.prop", "div.sqrt", "div.log"),
  div.index = c("richness", "shannon", "simpson", "mcintosh"),
  shannon.base = exp(1),
  div.fun = NULL,
  log.base = exp(1),
  metric = c("pooled", "mean"),
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

- qualitative:

  Name of columns with the qualitative traits as a character vector.

- method:

  The allocation method. Either `"div"` for constant or `"div.prop"` for
  proportional or `"div.log"` for logarithmic or `"div.sqrt"` for square
  root allocation.

- div.index:

  The diversity index to be used to estimate within cluster/group
  diversity.

- shannon.base:

  The logarithm base to be used for estimation of Shannon diversity
  index. Default is `exp(1)`.

- div.fun:

  A function to estimate diversity index from a factor vector of
  qualitative trait data.

- log.base:

  The logarithm base to be used for logarithmic method of sampling.
  Default is `exp(1)`.

- metric:

  The metric to be computed from the diversity index. Either `"pooled"`
  or `"mean"`.

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
in combination with the size of the cluster/group (See **`Methods`**).

There are several methods proposed on the basis of diversity indices
such as genetic multiplicity (G) dependent method based on the range of
genetic diversity (Yonezawa et al. 1995) , H strategy based on Nei's
gene diversity (Nei 1973) and a method based on the pooled Shannon
diversity index (Bisht et al. 1999; Mahajan et al. 1999) . Similarly,
measures such as expected proportion of heterozygous loci per individual
and effective number of alleles have also been employed as a diversity
measure for determining sample size (Franco et al. 2006) .

The within-cluster/group diversity is estimated as either pooled or mean
value of cluster/group-wise diversity indices. The following diversity
indices are implemented in this function.

- Shannon or Shannon-Weaver or Shannon-Wiener Diversity Index or Shannon
  entropy (\\H\\) (Shannon and Weaver 1949; Peet 1974)

- Simpson's Index of Diversity or Gini's Diversity Index or Gini-Simpson
  Index or Nei's Diversity Index or Nei's Variation Index (\\D\\) or
  Hurlbert’s probability of interspecific encounter (\\PIE\\) (Gini
  1912, 1912; Greenberg 1956; Berger and Parker 1970; Hurlbert 1971; Nei
  1973; Peet 1974)

- McIntosh Diversity Index (\\D\_{Mc}\\) (McIntosh 1967; Peet 1974)

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

## References

Berger WH, Parker FL (1970). “Diversity of planktonic foraminifera in
deep-sea sediments.” *Science*, **168**(3937), 1345–1347.  
  
Bisht IS, Mahajan RK, Gautam PL (1999). “Assessment of genetic
diversity, stratification of germplasm accessions in diversity groups
and sampling strategies for establishing a core collection of Indian
sesame (*Sesamum indicum* L.).” *Plant Genetic Resources Newsletter*,
**199 Supp.**, 35–46.  
  
Franco J, Crossa J, Warburton ML, Taba S (2006). “Sampling strategies
for conserving maize diversity when forming core subsets using genetic
markers.” *Crop Science*, **46**(2), 854–864.  
  
Gini C (1912). *Variabilita e Mutabilita. Contributo allo Studio delle
Distribuzioni e delle Relazioni Statistiche. \[Fasc. I.\]*. Tipogr. di
P. Cuppini, Bologna.  
  
Gini C (1912). “Variabilita e mutabilita.” In Pizetti E, Salvemini T
(eds.), *Memorie di Metodologica Statistica*. Liberia Eredi Virgilio
Veschi, Roma, Italy.  
  
Greenberg JH (1956). “The measurement of linguistic diversity.”
*Language*, **32**(1), 109.  
  
Hurlbert SH (1971). “The nonconcept of species diversity: A critique and
alternative parameters.” *Ecology*, **52**(4), 577–586.  
  
Mahajan RK, Bisht IS, Gautam PL (1999). “Sampling strategies for
developing Indian sesame core collection.” *Indian Journal of Plant
Genetic Resources*, **12**(01), 1–9.  
  
McIntosh RP (1967). “An index of diversity and the relation of certain
concepts to diversity.” *Ecology*, **48**(3), 392–404.  
  
Nei M (1973). “Analysis of gene diversity in subdivided populations.”
*Proceedings of the National Academy of Sciences*, **70**(12),
3321–3323.  
  
Peet RK (1974). “The measurement of species diversity.” *Annual Review
of Ecology and Systematics*, **5**(1), 285–307.  
  
Peet RK (1974). “The measurement of species diversity.” *Annual Review
of Ecology and Systematics*, **5**(1), 285–307.  
  
Shannon CE, Weaver W (1949). *The Mathematical Theory of Communication*,
number v. 2 in The Mathematical Theory of Communication. University of
Illinois Press.  
  
Yonezawa K, Nomura T, Morishima H (1995). “Sampling strategies for use
in stratified germplasm collections.” In Hodkin T, Brown ADH, van Hintum
TJL, Morales EAV (eds.), *Core Collections of Plant Genetic Resources*,
35–53. John Wiley & Sons, New York. ISBN 0-471-95545-0.

## See also

[`allocate.basic`](https://aravind-j.github.io/SampleCore/reference/allocate.basic.md),
[`allocate.distance`](https://aravind-j.github.io/SampleCore/reference/allocate.distance.md)

## Examples

``` r
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare example data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get data
data("cassava_EC_gp")
data <- cbind(genotypes = rownames(cassava_EC_gp), cassava_EC_gp)
row.names(data) <- NULL

# Column names of traits
quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW",
           "AVPW", "ARSR", "SRDM")
qual <- c("CUAL", "LNGS", "PTLC", "DSTA", "LFRT", "LBTEF", "CBTR", "NMLB",
          "ANGB", "CUAL9M", "LVC9M", "TNPR9M", "PL9M", "STRP", "STRC",
          "PSTR")

# Convert qualitative data columns to factor
data[, qual] <- lapply(data[, qual], as.factor)

# Convert quantitative data columns to qualitative scores
quant_to_score5 <- function(x) {

  brks <- unique( quantile(x,
                           probs = seq(0, 1, 0.2),
                           na.rm = TRUE))
  cut(x, breaks = brks,
      include.lowest = TRUE,
      labels = seq_len(length(brks) - 1))
}

data[, quant] <- lapply(data[, quant], quant_to_score5)

traits <- c(quant, qual)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Custom diversity index functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

div_fun_brillouin <- function(x) {
  n <- tabulate(x)
  n <- n[n > 0]
  N <- sum(n)
  if (N <= 1) {
    return(0)
  }
  (lgamma(N + 1) - sum(lgamma(n + 1)))/N
}

div_fun_margalef <- function(x) {
  tab <- tabulate(x)
  tab <- tab[tab > 0]
  S <- length(tab)
  N <- length(x)
  if (N <= 1) {
    return(0)
  }
  (S - 1)/log(N)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity allocation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shannon-Weaver Diversity Index
div_out_shannon1 <-
  allocate.diversity(data = data, names = "genotypes",
                    group = "Cluster",
                    qualitative = traits,
                    method = "div",
                    div.index = "shannon", metric = "pooled",
                    size = 0.2)
div_out_shannon1
#>   I  II III  IV   V  VI 
#>  63  53  49  59  55  59 

div_out_shannon2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.index = "shannon", metric = "mean",
                     size = 0.2)
div_out_shannon2
#>   I  II III  IV   V  VI 
#>  63  53  49  59  55  59 

##  Gini-Simpson Index
div_out_simpson1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.index = "simpson", metric = "pooled",
                     size = 0.2)
div_out_simpson1
#>   I  II III  IV   V  VI 
#>  62  53  49  59  55  58 

div_out_simpson2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.index = "simpson", metric = "mean",
                     size = 0.2)
div_out_simpson2
#>   I  II III  IV   V  VI 
#>  62  53  49  59  55  58 

## McIntosh Diversity Index
div_out_mcintosh1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.index = "mcintosh", metric = "pooled",
                     size = 0.2)
div_out_mcintosh1
#>   I  II III  IV   V  VI 
#>  63  54  49  59  54  59 

div_out_mcintosh2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.index = "mcintosh", metric = "mean",
                     size = 0.2)
div_out_mcintosh2
#>   I  II III  IV   V  VI 
#>  63  54  49  59  54  59 

## Richness
div_out_richness1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.index = "richness", metric = "pooled",
                     size = 0.2)
div_out_richness1
#>   I  II III  IV   V  VI 
#>  59  55  50  59  58  57 

div_out_richness2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.index = "richness", metric = "mean",
                     size = 0.2)
div_out_richness2
#>   I  II III  IV   V  VI 
#>  59  55  50  59  58  57 

## Brillouin Diversity Index
div_out_brillouin1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.fun = div_fun_brillouin, metric = "pooled",
                     size = 0.2)
div_out_brillouin1
#>   I  II III  IV   V  VI 
#>  63  53  48  59  55  58 

div_out_brillouin2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.fun = div_fun_brillouin, metric = "mean",
                     size = 0.2)
div_out_brillouin2
#>   I  II III  IV   V  VI 
#>  63  53  48  59  55  58 

## Margalef's richness Index
div_out_margalef1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.fun = div_fun_margalef, metric = "pooled",
                     size = 0.2)
div_out_margalef1
#>   I  II III  IV   V  VI 
#>  58  57  51  56  54  61 

div_out_margalef2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div",
                     div.fun = div_fun_margalef, metric = "mean",
                     size = 0.2)
div_out_margalef2
#>   I  II III  IV   V  VI 
#>  58  57  51  56  54  61 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity allocation & Proportional
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shannon-Weaver Diversity Index
dist_prop_out_shannon1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.index = "shannon", metric = "pooled",
                     size = 0.2)
dist_prop_out_shannon1
#>   I  II III  IV   V  VI 
#>  67  39  32  84  77  38 

dist_prop_out_shannon2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.index = "shannon", metric = "mean",
                     size = 0.2)
dist_prop_out_shannon2
#>   I  II III  IV   V  VI 
#>  67  39  32  84  77  38 

##  Gini-Simpson Index
dist_prop_out_simpson1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.index = "simpson", metric = "pooled",
                     size = 0.2)
dist_prop_out_simpson1
#>   I  II III  IV   V  VI 
#>  66  39  32  84  78  37 

dist_prop_out_simpson2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.index = "simpson", metric = "mean",
                     size = 0.2)
dist_prop_out_simpson2
#>   I  II III  IV   V  VI 
#>  66  39  32  84  78  37 

## McIntosh Diversity Index
dist_prop_out_mcintosh1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.index = "mcintosh", metric = "pooled",
                     size = 0.2)
dist_prop_out_mcintosh1
#>   I  II III  IV   V  VI 
#>  67  39  32  84  76  38 

dist_prop_out_mcintosh2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.index = "mcintosh", metric = "mean",
                     size = 0.2)
dist_prop_out_mcintosh2
#>   I  II III  IV   V  VI 
#>  67  39  32  84  76  38 

## Richness
div_out_richness1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "richness", metric = "pooled",
                     size = 0.2)
div_out_richness1
#>   I  II III  IV   V  VI 
#>  60  52  47  63  62  53 

div_out_richness2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "richness", metric = "mean",
                     size = 0.2)
div_out_richness2
#>   I  II III  IV   V  VI 
#>  60  52  47  63  62  53 

## Brillouin Diversity Index
div_out_brillouin1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.fun = div_fun_brillouin, metric = "pooled",
                     size = 0.2)
div_out_brillouin1
#>   I  II III  IV   V  VI 
#>  67  39  31  85  78  37 

div_out_brillouin2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.fun = div_fun_brillouin, metric = "mean",
                     size = 0.2)
div_out_brillouin2
#>   I  II III  IV   V  VI 
#>  67  39  31  85  78  37 

## Margalef's richness Index
div_out_margalef1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.fun = div_fun_margalef, metric = "pooled",
                     size = 0.2)
div_out_margalef1
#>   I  II III  IV   V  VI 
#>  63  42  34  81  77  40 

div_out_margalef2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.prop",
                     div.fun = div_fun_margalef, metric = "mean",
                     size = 0.2)
div_out_margalef2
#>   I  II III  IV   V  VI 
#>  63  42  34  81  77  40 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity allocation & Logarithmic
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shannon-Weaver Diversity Index
dist_log_out_shannon1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "shannon", metric = "pooled",
                     size = 0.2)
dist_log_out_shannon1
#>   I  II III  IV   V  VI 
#>  64  51  45  63  59  55 

dist_log_out_shannon2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "shannon", metric = "mean",
                     size = 0.2)
dist_log_out_shannon2
#>   I  II III  IV   V  VI 
#>  64  51  45  63  59  55 

##  Gini-Simpson Index
dist_log_out_simpson1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "simpson", metric = "pooled",
                     size = 0.2)
dist_log_out_simpson1
#>   I  II III  IV   V  VI 
#>  63  51  46  63  59  54 

dist_log_out_simpson2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "simpson", metric = "mean",
                     size = 0.2)
dist_log_out_simpson2
#>   I  II III  IV   V  VI 
#>  63  51  46  63  59  54 

## McIntosh Diversity Index
dist_log_out_mcintosh1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "mcintosh", metric = "pooled",
                     size = 0.2)
dist_log_out_mcintosh1
#>   I  II III  IV   V  VI 
#>  64  51  46  63  58  55 

dist_log_out_mcintosh2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "mcintosh", metric = "mean",
                     size = 0.2)
dist_log_out_mcintosh2
#>   I  II III  IV   V  VI 
#>  64  51  46  63  58  55 

## Richness
div_out_richness1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "richness", metric = "pooled",
                     size = 0.2)
div_out_richness1
#>   I  II III  IV   V  VI 
#>  60  52  47  63  62  53 

div_out_richness2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.index = "richness", metric = "mean",
                     size = 0.2)
div_out_richness2
#>   I  II III  IV   V  VI 
#>  60  52  47  63  62  53 

## Brillouin Diversity Index
div_out_brillouin1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.fun = div_fun_brillouin, metric = "pooled",
                     size = 0.2)
div_out_brillouin1
#>   I  II III  IV   V  VI 
#>  64  51  45  64  59  54 

div_out_brillouin2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.fun = div_fun_brillouin, metric = "mean",
                     size = 0.2)
div_out_brillouin2
#>   I  II III  IV   V  VI 
#>  64  51  45  64  59  54 

## Margalef's richness Index
div_out_margalef1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.fun = div_fun_margalef, metric = "pooled",
                     size = 0.2)
div_out_margalef1
#>   I  II III  IV   V  VI 
#>  59  54  48  60  58  57 

div_out_margalef2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.log",
                     div.fun = div_fun_margalef, metric = "mean",
                     size = 0.2)
div_out_margalef2
#>   I  II III  IV   V  VI 
#>  59  54  48  60  58  57 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity allocation & Square root
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shannon-Weaver Diversity Index
dist_sqrt_out_shannon1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.index = "shannon", metric = "pooled",
                     size = 0.2)
dist_sqrt_out_shannon1
#>   I  II III  IV   V  VI 
#>  66  46  40  71  66  48 

dist_sqrt_out_shannon2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.index = "shannon", metric = "mean",
                     size = 0.2)
dist_sqrt_out_shannon2
#>   I  II III  IV   V  VI 
#>  66  46  40  71  66  48 

##  Gini-Simpson Index
dist_sqrt_out_simpson1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.index = "simpson", metric = "pooled",
                     size = 0.2)
dist_sqrt_out_simpson1
#>   I  II III  IV   V  VI 
#>  65  46  40  71  67  47 

dist_sqrt_out_simpson2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.index = "simpson", metric = "mean",
                     size = 0.2)
dist_sqrt_out_simpson2
#>   I  II III  IV   V  VI 
#>  65  46  40  71  67  47 

## McIntosh Diversity Index
dist_sqrt_out_mcintosh1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.index = "mcintosh", metric = "pooled",
                     size = 0.2)
dist_sqrt_out_mcintosh1
#>   I  II III  IV   V  VI 
#>  66  47  40  71  65  48 

dist_sqrt_out_mcintosh2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.index = "mcintosh", metric = "mean",
                     size = 0.2)
dist_sqrt_out_mcintosh2
#>   I  II III  IV   V  VI 
#>  66  47  40  71  65  48 

## Richness
div_out_richness1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.index = "richness", metric = "pooled",
                     size = 0.2)
div_out_richness1
#>   I  II III  IV   V  VI 
#>  61  47  41  72  69  46 

div_out_richness2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.index = "richness", metric = "mean",
                     size = 0.2)
div_out_richness2
#>   I  II III  IV   V  VI 
#>  61  47  41  72  69  46 

## Brillouin Diversity Index
div_out_brillouin1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.fun = div_fun_brillouin, metric = "pooled",
                     size = 0.2)
div_out_brillouin1
#>   I  II III  IV   V  VI 
#>  66  46  39  72  67  47 

div_out_brillouin2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.fun = div_fun_brillouin, metric = "mean",
                     size = 0.2)
div_out_brillouin2
#>   I  II III  IV   V  VI 
#>  66  46  39  72  67  47 

## Margalef's richness Index
div_out_margalef1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.fun = div_fun_margalef, metric = "pooled",
                     size = 0.2)
div_out_margalef1
#>   I  II III  IV   V  VI 
#>  61  49  42  68  65  50 

div_out_margalef2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = traits,
                     method = "div.sqrt",
                     div.fun = div_fun_margalef, metric = "mean",
                     size = 0.2)
div_out_margalef2
#>   I  II III  IV   V  VI 
#>  61  49  42  68  65  50 
```
