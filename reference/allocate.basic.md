# Allocation of Entries to be Selected from Clusters/Groups based on Size for Core Collection Development

Estimate the number of entries to be allocated from each cluster/group
in the entire collection to construct a core collection on the basis of
cluster/group size. The following strategies are implemented.

- Constant

- Proportional

- Logarithmic

- Square root

The different methods to determine the number of entries from each group
or clusters implemented in `allocate.basic` are as follows.

## Usage

``` r
allocate.basic(
  data,
  names,
  group,
  method = c("const", "prop", "log", "sqrt"),
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

- method:

  The allocation method. Either `"const"` for constant or `"prop"` for
  proportional or `"log"` for logarithmic or `"sqrt"` for square root
  allocation.

- log.base:

  The logarithm base to be used for logarithmic method of sampling.
  Default is `exp(1)`.

- size:

  The desired core set size proportion.

## Details

These are different methods which estimate the number of entries only on
the basis of total number of entries in each cluster/group.

Brown (1989) proposed the constant (C), proportional (P) and logarithmic
(L) methods and later a similar square root method was proposed by
Huaman et al. (1999) .

### Constant method

From an entire collection of size \\N\\, to construct a core set of
sample size \\n\\, the number of entries to be selected from the \\i\\th
group among \\1 \cdots g\\ groups (\\n\_{i}\\) is estimated as below.

\\n\_{i} = \frac{n}{g} \times N\\

### Proportional method

Here the number of entries to be selected is proportional to the
cluster/group size (\\N\_{i}\\) as below.

\\n\_{i} = n \times \frac{N\_{i}}{\sum\_{i=1}^{g}N\_{i}}\\

\\n\_{i} = n \times \frac{N\_{i}}{N}\\

### Logarithmic method

Here the number of entries to be selected is proportional to the
logarithm of the cluster/group size (\\N\_{i}\\) as below.

\\n\_{i} = n \times
\frac{\log{(N\_{i})}}{\sum\_{i=1}^{g}\log{(N\_{i})}}\\

### Square root method

Here the number of entries to be selected is proportional to the square
root of the cluster/group size (\\N\_{i}\\) as below.

\\n\_{i} = n \times \frac{\sqrt{N\_{i}}}{\sum\_{i=1}^{g}\sqrt{N\_{i}}}\\

## References

Brown AHD (1989). “Core collections: A practical approach to genetic
resources management.” *Genome*, **31**(2), 818–824.  
  
Huaman Z, Aguilar C, Ortiz R (1999). “Selecting a Peruvian sweetpotato
core collection on the basis of morphological, eco-geographical, and
disease and pest reaction data:.” *Theoretical and Applied Genetics*,
**98**(5), 840–844.
