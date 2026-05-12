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
  quantitative,
  qualitative,
  method = c("div", "div.prop", "div.sqrt", "div.log"),
  div.index = c("shannon", "simpson", "mcintosh"),
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

- quantitative:

  Name of columns with the quantitative traits as a character vector.

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
  
Hurlbert SH (1971). “The nonconcept of species diversity: a critique and
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

## Examples

``` r
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare example data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity allocation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shannon-Weaver Diversity Index
dist_out_shannon1 <-
  allocate.diversity(data = data, names = "genotypes",
                    group = "Cluster",
                    qualitative = qual, quantitative = quant,
                    method = "div",
                    div.index = "shannon", metric = "pooled",
                    size = 0.2)
dist_out_shannon1
#>   I  II III  IV   V  VI 
#>  67  51  46  60  54  60 

dist_out_shannon2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div",
                     div.index = "shannon", metric = "mean",
                     size = 0.2)
dist_out_shannon2
#>   I  II III  IV   V  VI 
#>  67  51  46  60  54  60 

##  Gini-Simpson Index
dist_out_simpson1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div",
                     div.index = "simpson", metric = "pooled",
                     size = 0.2)
dist_out_simpson1
#>   I  II III  IV   V  VI 
#>  66  51  46  60  55  59 

dist_out_simpson2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div",
                     div.index = "simpson", metric = "mean",
                     size = 0.2)
dist_out_simpson2
#>   I  II III  IV   V  VI 
#>  66  51  46  60  55  59 

## McIntosh Diversity Index
dist_out_mcintosh1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div",
                     div.index = "mcintosh", metric = "pooled",
                     size = 0.2)
dist_out_mcintosh1
#>   I  II III  IV   V  VI 
#>  66  52  47  60  53  59 

dist_out_mcintosh2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div",
                     div.index = "mcintosh", metric = "mean",
                     size = 0.2)
dist_out_mcintosh2
#>   I  II III  IV   V  VI 
#>  66  52  47  60  53  59 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity allocation & Proportional
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shannon-Weaver Diversity Index
dist_prop_out_shannon1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.prop",
                     div.index = "shannon", metric = "pooled",
                     size = 0.2)
dist_prop_out_shannon1
#>   I  II III  IV   V  VI 
#>  71  37  30  85  75  39 

dist_prop_out_shannon2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.prop",
                     div.index = "shannon", metric = "mean",
                     size = 0.2)
dist_prop_out_shannon2
#>   I  II III  IV   V  VI 
#>  71  37  30  85  75  39 

##  Gini-Simpson Index
dist_prop_out_simpson1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.prop",
                     div.index = "simpson", metric = "pooled",
                     size = 0.2)
dist_prop_out_simpson1
#>   I  II III  IV   V  VI 
#>  70  37  30  86  77  38 

dist_prop_out_simpson2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.prop",
                     div.index = "simpson", metric = "mean",
                     size = 0.2)
dist_prop_out_simpson2
#>   I  II III  IV   V  VI 
#>  70  37  30  86  77  38 

## McIntosh Diversity Index
dist_prop_out_mcintosh1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.prop",
                     div.index = "mcintosh", metric = "pooled",
                     size = 0.2)
dist_prop_out_mcintosh1
#>   I  II III  IV   V  VI 
#>  71  38  30  86  74  38 

dist_prop_out_mcintosh2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.prop",
                     div.index = "mcintosh", metric = "mean",
                     size = 0.2)
dist_prop_out_mcintosh2
#>   I  II III  IV   V  VI 
#>  71  38  30  86  74  38 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity allocation & Logarithmic
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shannon-Weaver Diversity Index
dist_log_out_shannon1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.log",
                     div.index = "shannon", metric = "pooled",
                     size = 0.2)
dist_log_out_shannon1
#>   I  II III  IV   V  VI 
#>  68  49  42  64  58  56 

dist_log_out_shannon2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.log",
                     div.index = "shannon", metric = "mean",
                     size = 0.2)
dist_log_out_shannon2
#>   I  II III  IV   V  VI 
#>  68  49  42  64  58  56 

##  Gini-Simpson Index
dist_log_out_simpson1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.log",
                     div.index = "simpson", metric = "pooled",
                     size = 0.2)
dist_log_out_simpson1
#>   I  II III  IV   V  VI 
#>  67  49  43  64  59  55 

dist_log_out_simpson2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.log",
                     div.index = "simpson", metric = "mean",
                     size = 0.2)
dist_log_out_simpson2
#>   I  II III  IV   V  VI 
#>  67  49  43  64  59  55 

## McIntosh Diversity Index
dist_log_out_mcintosh1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.log",
                     div.index = "mcintosh", metric = "pooled",
                     size = 0.2)
dist_log_out_mcintosh1
#>   I  II III  IV   V  VI 
#>  68  49  44  65  57  55 

dist_log_out_mcintosh2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.log",
                     div.index = "mcintosh", metric = "mean",
                     size = 0.2)
dist_log_out_mcintosh2
#>   I  II III  IV   V  VI 
#>  68  49  44  65  57  55 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Diversity allocation & Square root
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Shannon-Weaver Diversity Index
dist_sqrt_out_shannon1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.sqrt",
                     div.index = "shannon", metric = "pooled",
                     size = 0.2)
dist_sqrt_out_shannon1
#>   I  II III  IV   V  VI 
#>  70  44  37  72  65  49 

dist_sqrt_out_shannon2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.sqrt",
                     div.index = "shannon", metric = "mean",
                     size = 0.2)
dist_sqrt_out_shannon2
#>   I  II III  IV   V  VI 
#>  70  44  37  72  65  49 

##  Gini-Simpson Index
dist_sqrt_out_simpson1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.sqrt",
                     div.index = "simpson", metric = "pooled",
                     size = 0.2)
dist_sqrt_out_simpson1
#>   I  II III  IV   V  VI 
#>  69  44  38  73  66  48 

dist_sqrt_out_simpson2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.sqrt",
                     div.index = "simpson", metric = "mean",
                     size = 0.2)
dist_sqrt_out_simpson2
#>   I  II III  IV   V  VI 
#>  69  44  38  73  66  48 

## McIntosh Diversity Index
dist_sqrt_out_mcintosh1 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.sqrt",
                     div.index = "mcintosh", metric = "pooled",
                     size = 0.2)
dist_sqrt_out_mcintosh1
#>   I  II III  IV   V  VI 
#>  69  45  38  73  64  48 

dist_sqrt_out_mcintosh2 <-
  allocate.diversity(data = data, names = "genotypes",
                     group = "Cluster",
                     qualitative = qual, quantitative = quant,
                     method = "div.sqrt",
                     div.index = "mcintosh", metric = "mean",
                     size = 0.2)
dist_sqrt_out_mcintosh2
#>   I  II III  IV   V  VI 
#>  69  45  38  73  64  48 
```
