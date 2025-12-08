# Diversity Indices

Compute the following diversity indices.

- Simpson's and related indices

  - Simpson's Index (\\d\\) (Simpson 1949; Peet 1974)

  - Simpson's Index of Diversity or Gini's Diversity Index or
    Gini-Simpson Index or Nei's Diversity Index or Nei's Variation Index
    (\\D\\) (Gini 1912, 1912; Greenberg 1956; Berger and Parker 1970;
    Nei 1973; Peet 1974)

  - Maximum Simpson's Index of Diversity or Maximum Nei's
    Diversity/Variation Index (\\D\_{max}\\) (Hennink and Zeven 1990)

  - Simpson's Reciprocal Index or Hill's \\N\_{2}\\ (\\D\_{R}\\)
    (Williams 1964; Hill 1973)

  - Relative Simpson's Index of Diversity or Relative Nei's
    Diversity/Variation Index (\\D'\\) (Hennink and Zeven 1990)

- Shannon-Weaver and related indices

  - Shannon or Shannon-Weaver or Shannon-Weiner Diversity Index (\\H\\)
    (Shannon and Weaver 1949; Peet 1974)

  - Maximum Shannon-Weaver Diversity Index (\\H\_{max}\\) (Hennink and
    Zeven 1990)

  - Relative Shannon-Weaver Diversity Index or Shannon Equitability
    Index (\\H'\\) (Hennink and Zeven 1990)

- McIntosh Diversity Index

  - McIntosh Diversity Index (\\D\_{Mc}\\) (McIntosh 1967; Peet 1974)

## Usage

``` r
diversity.indices(
  x,
  index = c("simpson", "simpson.cmp", "simpson.max", "simpson.inv", "simpson.rel",
    "shannon", "shannon.max", "shannon.rel", "mcintosh"),
  base = 2
)
```

## Arguments

- x:

  The qualitative trait data as a vector of class factor.

- base:

  The logarithm base to be used for computation of Shannon-Weaver
  Diversity Index (\\I\\). Default is 2.

## Value

The diversity index value.

## Details

The diversity indices and the corresponding statistical tests
implemented in `diversity.indices` are as follows.

### Simpson's and related indices

Simpson's index (\\d\\) which estimates the probability that two
accessions randomly selected will belong to the same class of a
qualitative trait, is computed as follows (Simpson 1949; Peet 1974) .

\\d = \sum\_{i = 1}^{k}p\_{i}^{2}\\

Where, \\p\_{i}\\ denotes the proportion/fraction/frequency of
accessions in the \\i\\th class of a qualitative trait and \\k\\ is the
number of classes for the qualitative trait.

The value of \\d\\ can range from 0 to 1 with 0 representing maximum
diversity and 1, no diversity.

The complement of \\d\\ (\\d\\ is subtracted from 1) is called the
Simpson's index of diversity (\\D\\) (Greenberg 1956; Berger and Parker
1970; Peet 1974; Hennink and Zeven 1990) originally suggested by Gini
(1912, 1912) and described in literature as Gini's diversity index or
Gini-Simpson index. It is the same as Nei's diversity index or Nei's
variation index (Nei 1973; Hennink and Zeven 1990) . Greater the value
of \\D\\, greater the diversity with a range from 0 to 1.

\\D = 1 - d\\

The maximum value of \\D\\, \\D\_{max}\\ occurs when accessions are
uniformly distributed across the classes in the qualitative trait and is
computed as follows (Hennink and Zeven 1990) .

\\D\_{max} = 1 - \frac{1}{k}\\

Reciprocal of \\d\\ gives the Simpson's reciprocal index (\\D\_{R}\\)
(Williams 1964; Hennink and Zeven 1990) and can range from 1 to \\k\\.
This was also described in Hill (1973) as (\\N\_{2}\\).

\\D\_{R} = \frac{1}{d}\\

Relative Simpson's index of diversity or Relative Nei's
diversity/variation index (\\H'\\) (Hennink and Zeven 1990) is defined
as follows (Peet 1974) .

\\D' = \frac{D}{D\_{max}}\\

### Shannon-Weaver and related indices

An index of information \\H\\, was described by Shannon and Weaver
(1949) as follows.

\\H = -\sum\_{i=1}^{k}p\_{i} \log\_{2}(p\_{i})\\

\\H\\ is described as Shannon or Shannon-Weaver or Shannon-Weiner
diversity index in literature.

Alternatively, \\H\\ is also computed using natural logarithm instead of
logarithm to base 2.

\\H = -\sum\_{i=1}^{k}p\_{i} \ln(p\_{i})\\

The maximum value of \\H\\ (\\H\_{max}\\) is \\\ln(k)\\. This value
occurs when each class for a qualitative trait has the same proportion
of accessions.

\\H\_{max} = \log\_{2}(k)\\\\ \textrm{OR} \\\\ H\_{max} = \ln(k)\\

The relative Shannon-Weaver diversity index or Shannon equitability
index (\\H'\\) is the Shannon diversity index (\\I\\) divided by the
maximum diversity (\\H\_{max}\\).

\\H' = \frac{H}{H\_{max}}\\

### McIntosh Diversity Index

A similar index of diversity was described by McIntosh (1967) as follows
(\\D\_{Mc}\\) (Peet 1974) .

\\D\_{Mc} = \frac{N - \sqrt{\sum\_{i=1}^{k}n\_{i}^2}}{N - \sqrt{N}}\\

Where, \\n\_{i}\\ denotes the number of accessions in the \\i\\th class
for a qualitative trait and \\N\\ is the total number of accessions so
that \\p\_{i} = {n\_{i}}/{N}\\.

## References

Berger WH, Parker FL (1970). “Diversity of planktonic foraminifera in
deep-sea sediments.” *Science*, **168**(3937), 1345–1347.  
  
Gini C (1912). *Variabilita e Mutabilita. Contributo allo Studio delle
Distribuzioni e delle Relazioni Statistiche. \[Fasc. I.\]*. Tipogr. di
P. Cuppini, Bologna.  
  
Gini C (1912). “Variabilita e mutabilita.” In Pizetti E, Salvemini T
(eds.), *Memorie di Metodologica Statistica*. Liberia Eredi Virgilio
Veschi, Roma, Italy.  
  
Greenberg JH (1956). “The measurement of linguistic diversity.”
*Language*, **32**(1), 109.  
  
Hennink S, Zeven AC (1990). “The interpretation of Nei and
Shannon-Weaver within population variation indices.” *Euphytica*,
**51**(3), 235–240.  
  
Hill MO (1973). “Diversity and evenness: A unifying notation and its
consequences.” *Ecology*, **54**(2), 427–432.  
  
McIntosh RP (1967). “An index of diversity and the relation of certain
concepts to diversity.” *Ecology*, **48**(3), 392–404.  
  
Nei M (1973). “Analysis of gene diversity in subdivided populations.”
*Proceedings of the National Academy of Sciences*, **70**(12),
3321–3323.  
  
Peet RK (1974). “The measurement of species diversity.” *Annual Review
of Ecology and Systematics*, **5**(1), 285–307.  
  
Shannon CE, Weaver W (1949). *The Mathematical Theory of Communication*,
number v. 2 in The Mathematical Theory of Communication. University of
Illinois Press.  
  
Simpson EH (1949). “Measurement of diversity.” *Nature*, **163**(4148),
688–688.  
  
Williams CB (1964). *Patterns in the Balance of Nature and Related
Problems in Quantitative Ecology*. Academic Press.
