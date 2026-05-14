# Selection of Entries from Clusters/Groups by Random Sampling

Select entries from cluster/groups in the entire collection by random
sampling according to allocation specified.

## Usage

``` r
select.random(data, names, group, alloc, always.selected = NULL)
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

- always.selected:

  Names of accessions to be always included in the core set as a
  character vector.

## Value

A named list where each element contains the selected entry identifiers
for a cluster/group.

## Details

For each cluster/group entries are selected randomly according to the
allocation provided (Brown 1989; Brown and van Hintum 2000) . Entries
listed as `always.selected` are mandatorily included in the selection.
Warnings are issued if requested allocation is smaller than the number
of always-selected entries in a cluster/group and/or when the
cluster/group does not contain enough remaining entries to fulfill the
allocation.

## References

Brown AHD (1989). “Core collections: A practical approach to genetic
resources management.” *Genome*, **31**(2), 818–824.  
  
Brown AHD, van Hintum TJL (2000). *Core Collections of Plant Genetic
Resources*. Bioversity International. ISBN 92-9043-454-6.

## Examples

``` r

library(cluster)

# Get data
data(cassava_EC_gp)

data <- cbind(genotypes = rownames(cassava_EC_gp), cassava_EC_gp)
row.names(data) <- NULL

# Prepare inputs
counts <- c(I = 31, II = 31, III = 18, IV = 35, V = 40, VI = 17)

mand_accns <-
  c("TMe-34", "TMe-3423", "TMe-2018", "TMe-801", "TMe-551")

# Specify the seed
set.seed(123)

# Fetch selected accessions
sel_random_out <-
  select.random(data = data, names = "genotypes",
                group = "Cluster", alloc = counts,
                always.selected = mand_accns)

sel_random_out
#> $I
#>  [1] "TMe-3423" "TMe-2906" "TMe-489"  "TMe-2984" "TMe-3262" "TMe-3314"
#>  [7] "TMe-3418" "TMe-3667" "TMe-896"  "TMe-2237" "TMe-2253" "TMe-3471"
#> [13] "TMe-2993" "TMe-3396" "TMe-486"  "TMe-815"  "TMe-299"  "TMe-3065"
#> [19] "TMe-3460" "TMe-3149" "TMe-2027" "TMe-2066" "TMe-1344" "TMe-694" 
#> [25] "TMe-990"  "TMe-3110" "TMe-3514" "TMe-756"  "TMe-264"  "TMe-3266"
#> [31] "TMe-1581"
#> 
#> $II
#>  [1] "TMe-1172" "TMe-2077" "TMe-2127" "TMe-2204" "TMe-2021" "TMe-1461"
#>  [7] "TMe-3009" "TMe-2891" "TMe-2414" "TMe-1242" "TMe-601"  "TMe-1323"
#> [13] "TMe-3258" "TMe-2866" "TMe-1919" "TMe-477"  "TMe-40"   "TMe-2978"
#> [19] "TMe-289"  "TMe-2797" "TMe-3531" "TMe-1250" "TMe-2903" "TMe-53"  
#> [25] "TMe-1623" "TMe-960"  "TMe-85"   "TMe-433"  "TMe-3805" "TMe-431" 
#> [31] "TMe-3800"
#> 
#> $III
#>  [1] "TMe-1910" "TMe-425"  "TMe-3071" "TMe-3326" "TMe-785"  "TMe-1443"
#>  [7] "TMe-3544" "TMe-3750" "TMe-1176" "TMe-2901" "TMe-3324" "TMe-174" 
#> [13] "TMe-3679" "TMe-3321" "TMe-2270" "TMe-2977" "TMe-2362" "TMe-3174"
#> 
#> $IV
#>  [1] "TMe-34"   "TMe-801"  "TMe-3615" "TMe-3204" "TMe-1297" "TMe-3191"
#>  [7] "TMe-2788" "TMe-594"  "TMe-1118" "TMe-3259" "TMe-3451" "TMe-317" 
#> [13] "TMe-1430" "TMe-3108" "TMe-87"   "TMe-2059" "TMe-81"   "TMe-3248"
#> [19] "TMe-1148" "TMe-2332" "TMe-2399" "TMe-279"  "TMe-3527" "TMe-608" 
#> [25] "TMe-3312" "TMe-3542" "TMe-1325" "TMe-1867" "TMe-1700" "TMe-3004"
#> [31] "TMe-3054" "TMe-1489" "TMe-735"  "TMe-3089" "TMe-3055"
#> 
#> $V
#>  [1] "TMe-2018" "TMe-551"  "TMe-167"  "TMe-1285" "TMe-800"  "TMe-832" 
#>  [7] "TMe-877"  "TMe-823"  "TMe-547"  "TMe-1300" "TMe-481"  "TMe-685" 
#> [13] "TMe-1332" "TMe-439"  "TMe-1101" "TMe-1715" "TMe-326"  "TMe-475" 
#> [19] "TMe-797"  "TMe-884"  "TMe-2530" "TMe-2139" "TMe-3332" "TMe-2820"
#> [25] "TMe-712"  "TMe-2121" "TMe-866"  "TMe-1037" "TMe-977"  "TMe-1453"
#> [31] "TMe-3311" "TMe-373"  "TMe-1733" "TMe-2790" "TMe-1522" "TMe-2055"
#> [37] "TMe-510"  "TMe-378"  "TMe-2319" "TMe-456" 
#> 
#> $VI
#>  [1] "TMe-1539" "TMe-2196" "TMe-514"  "TMe-907"  "TMe-1216" "TMe-1392"
#>  [7] "TMe-2389" "TMe-696"  "TMe-690"  "TMe-1995" "TMe-2535" "TMe-1676"
#> [13] "TMe-1110" "TMe-1302" "TMe-631"  "TMe-1017" "TMe-1007"
#> 

# Get distance matrix - Only for visualization
quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW",
           "AVPW", "ARSR", "SRDM")
qual <- c("CUAL", "LNGS", "PTLC", "DSTA", "LFRT", "LBTEF", "CBTR", "NMLB",
          "ANGB", "CUAL9M", "LVC9M", "TNPR9M", "PL9M", "STRP", "STRC",
          "PSTR")

# Convert qualitative data columns to factor
cassava_EC_gp[, qual] <- lapply(cassava_EC_gp[, qual], as.factor)

# Standardise quantitative data column
cassava_EC_gp[, quant] <- lapply(cassava_EC_gp[, quant], function(x) {
  scale(x)[, 1]
})

gp_vec <- setNames(as.character(data[, "Cluster"]), data[, "genotypes"])

# Get the Gower's distance matrix
dist_matrix <- daisy(x = cassava_EC_gp[, c(qual, quant)],
                     metric = "gower")

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(sel_random_out, use.names = FALSE))
```
