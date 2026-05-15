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

## Value

A named numeric vector specifying the number of entries to be selected
from each cluster/group. The vector names correspond to the levels of
the "`"group"` column, and values indicate the number of elements to be
selected from each level.

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

## See also

[`allocate.distance`](https://aravind-j.github.io/SampleCore/reference/allocate.distance.md),
[`allocate.diversity`](https://aravind-j.github.io/SampleCore/reference/allocate.diversity.md)

## Examples

``` r
# Get data
data("cassava_EC_gp")

set.seed(123)
cassava_EC_gp <- cassava_EC_gp[sample(1:nrow(cassava_EC_gp), 500), ]

data <- cassava_EC_gp

data <- cbind(genotypes = rownames(cassava_EC_gp), cassava_EC_gp)
row.names(data) <- NULL

# Constant allocation
const_out <-
  allocate.basic(data = data, names = "genotypes",
                 group = "Cluster", method = "const",
                 size = 0.2)
const_out
#>   I  II III  IV   V  VI 
#>  17  17  17  17  17  17 

# Proportional allocation
prop_out <-
  allocate.basic(data = data, names = "genotypes",
                 group = "Cluster", method = "prop",
                 size = 0.2)
prop_out
#>   I  II III  IV   V  VI 
#>  17  11  11  27  23   9 

# Logarithmic allocation
log_out <-
  allocate.basic(data = data, names = "genotypes",
                 group = "Cluster", method = "log",
                 size = 0.2)
log_out
#>   I  II III  IV   V  VI 
#>  17  16  15  19  18  15 

# Square root allocation
sqrt_out <-
  allocate.basic(data = data, names = "genotypes",
                 group = "Cluster", method = "sqrt",
                 size = 0.2)
sqrt_out
#>   I  II III  IV   V  VI 
#>  17  14  14  22  20  13 
```
