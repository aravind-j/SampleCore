# Sample Size Estimation

Estimate the number of accessions to be sampled from each group/cluster
in the entire collection to construct a core collection. The following
strategies are implemented.

- Basic methods (Brown 1989; Huaman et al. 1999)

  - Constant

  - Proportional

  - Logarithmic

  - Square root

- Diversity dependent methods (Yonezawa et al. 1995; Schoen and Brown
  1993; Bisht et al. 1999; Mahajan et al. 1999; Franco et al. 2005)

  - Diversity

  - Diversity & Proportional

  - Diversity & Logarithmic

  - Diversity & Square root

## Usage

``` r
sample.size(
  data,
  names,
  group,
  dist,
  quantitative,
  qualitative,
  size = 0.2,
  sample.method = c("const", "prop", "log", "sqrt", "diversity", "diverstiy.prop",
    "diversity.sqrt", "diversity.log")
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

- dist:

  A precomputed distance matrix of distance measures between the
  accessions in `data`.

- size:

  The desired core set size proportion.

- sample.method:

## Details

The different methods to determine the number of entries from each group
or clusters implemented in `sample.size` are as follows.

### Basic methods

These are different methods which estimate the number of entries only on
the basis of total number of accessions in each group/cluster. Brown
(1989) proposed the constant (C), proportional (P) and logarithmic (L)
methods and later a similar square root method was proposed by Huaman et
al. (1999) .

#### Constant method

From an entire collection of size \\N\\, to construct a core set of
sample size \\n\\, the number of entries to be selected from the \\i\\th
group among \\1 \cdots g\\ groups (\\n\_{i}\\) is estimated as below.

\\n\_{i} = \frac{n}{g} \times N\\

#### Proportional method

Here the number of entries to be selected is proportional to the
group/cluster size (\\N\_{i}\\) as below.

\\n\_{i} = n \times \frac{N\_{i}}{\sum\_{i=1}^{g}N\_{i}}\\

\\n\_{i} = n \times \frac{N\_{i}}{N}\\

#### Logarithmic method

Here the number of entries to be selected is proportional to the
logarithm of the group/cluster size (\\N\_{i}\\) as below.

\\n\_{i} = n \times
\frac{\log{(N\_{i})}}{\sum\_{i=1}^{g}\log{(N\_{i})}}\\

#### Square root method

Here the number of entries to be selected is proportional to the square
root of the group/cluster size (\\N\_{i}\\) as below.

\\n\_{i} = n \times \frac{\sqrt{N\_{i}}}{\sum\_{i=1}^{g}\sqrt{N\_{i}}}\\

### Diversity dependent methods

These are different methods which estimate the number of entries on the
basis of how diverse are the accessions in each group/cluster. There are
several methods proposed on the basis of diversity indices such as
genetic multiplicity (G) dependent method based on the range of genetic
diversity (Yonezawa et al. 1995) , H strategy based on Nei's gene
diversity (Nei 1973) and a method based on the pooled Shannon diversity
index (Bisht et al. 1999; Mahajan et al. 1999) . Similarly, measures
such as expected proportion of heterozygous loci per individual and
effective number of alleles have also been employed as a diversity
measure for determining sample size (Franco et al. 2006) . Franco et al.
(2005) proposed a method based on mean Gower's distance (Gower 1971)
which was also extended to other distance measure averages named D
Allocation strategy (Franco et al. 2006) . These methods were also
combined with the proportional and logarithmic methods. For example, the
GP and GL strategy of Bisht et al. (1999) and Mahajan et al. (1999) as
well as the NY and LD allocation methods of Franco et al. (2005) .

#### Diversity method

From an entire collection of size \\N\\, to construct a core set of
sample size \\n\\, the number of entries to be selected from the \\i\\th
group among \\1 \cdots g\\ groups (\\n\_{i}\\) is estimated as below.

\\n\_{i} = n \times \frac{D\_{i}}{\sum\_{i=1}^{g}D\_{i}}\\

Where, \\D\_{i}\\ is a measure of the extent of diversity present in the
\\i\\th cluster.

\\D\\ can be either 1) Range of a diversity index 2) Pooled value of a
diversity index or 3) Mean genetic distance.

#### Diversity and proportional method

Here the number of entries to be selected is proportional to the
diversity of the group/cluster (\\D\_{i}\\) weighted by the the
group/cluster size (\\N\_{i}\\).

\\n\_{i} = n \times \frac{N\_{i}D\_{i}}{\sum\_{i=1}^{g}N\_{i}D\_{i}}\\

#### Diversity and logarithmic method

Here the number of entries to be selected is proportional to the
diversity of the group/cluster (\\D\_{i}\\) weighted by the logarithm of
the group/cluster size (\\N\_{i}\\).

\\n\_{i} = n \times
\frac{\log(N\_{i})D\_{i}}{\sum\_{i=1}^{g}\log(N\_{i})D\_{i}}\\

#### Diversity and square root method

Here the number of entries to be selected is proportional to the
diversity of the group/cluster (\\D\_{i}\\) weighted by the square root
of the group/cluster size (\\N\_{i}\\).

\\n\_{i} = n \times
\frac{\sqrt{N\_{i}}D\_{i}}{\sum\_{i=1}^{g}\sqrt{N\_{i}}D\_{i}}\\

## References

Bisht IS, Mahajan RK, Gautam PL (1999). “Assessment of genetic
diversity, stratification of germplasm accessions in diversity groups
and sampling strategies for establishing a core collection of Indian
sesame (*Sesamum indicum* L.).” *Plant Genetic Resources Newsletter*,
**199 Supp.**, 35–46.  
  
Brown AHD (1989). “Core collections: A practical approach to genetic
resources management.” *Genome*, **31**(2), 818–824.  
  
Franco J, Crossa J, Taba S, Shands H (2005). “A sampling strategy for
conserving genetic diversity when forming core subsets.” *Crop Science*,
**45**(3), 1035–1044.  
  
Franco J, Crossa J, Warburton ML, Taba S (2006). “Sampling strategies
for conserving maize diversity when forming core subsets using genetic
markers.” *Crop Science*, **46**(2), 854–864.  
  
Gower JC (1971). “A general coefficient of similarity and some of its
properties.” *Biometrics*, **27**(4), 857–871.  
  
Huaman Z, Aguilar C, Ortiz R (1999). “Selecting a Peruvian sweetpotato
core collection on the basis of morphological, eco-geographical, and
disease and pest reaction data:.” *Theoretical and Applied Genetics*,
**98**(5), 840–844.  
  
Mahajan RK, Bisht IS, Gautam PL (1999). “Sampling strategies for
developing Indian sesame core collection.” *Indian Journal of Plant
Genetic Resources*, **12**(01), 1–9.  
  
Nei M (1973). “Analysis of gene diversity in subdivided populations.”
*Proceedings of the National Academy of Sciences*, **70**(12),
3321–3323.  
  
Schoen DJ, Brown AH (1993). “Conservation of allelic richness in wild
crop relatives is aided by assessment of genetic markers.” *Proceedings
of the National Academy of Sciences*, **90**(22), 10623–10627.  
  
Yonezawa K, Nomura T, Morishima H (1995). “Sampling strategies for use
in stratified germplasm collections.” In Hodkin T, Brown ADH, Hintum
TJLv, Morales EAV (eds.), *Core Collections of Plant Genetic Resources*,
35–53. John Wiley \\ Sons, New York. ISBN 0-471-95545-0.
