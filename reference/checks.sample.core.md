# Common checks for all functions

Not exported. Strictly internal

## Usage

``` r
checks.sample.core(
  data,
  names,
  size,
  group,
  quantitative = NULL,
  qualitative = NULL,
  dist.mat = NULL,
  log.base = NULL,
  always.selected = NULL,
  mode = C("alloc", "sel")
)
```

## Arguments

- data:

  The data as a data frame object. The data frame should possess one row
  per individual and columns with the individual names and multiple
  trait/character data.

- names:

  Name of column with the accession names as a character string.

- size:

  The desired core set size proportion.

- group:

  Name of column with the accession group/cluster names as a character
  string.

- quantitative:

  Name of columns with the quantitative traits as a character vector.

- qualitative:

  Name of columns with the qualitative traits as a character vector.

- dist.mat:

  A precomputed distance matrix of distance measures between the
  accessions in `data`.

- log.base:

  The logarithm base to be used for logarithmic method of sampling.
  Default is `exp(1)`.
