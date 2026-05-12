# Selection of Entries from Clusters/Groups by Random Sampling

Select entries from cluster/groups in the entire collection by random
sampling according to allocation specified.

## Usage

``` r
select.random(data, names, group, alloc, always.selected)
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
counts <- c(I = 61, II = 41, III = 37, IV = 81, V = 80, VI = 37)

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
#> [31] "TMe-1581" "TMe-3211" "TMe-1940" "TMe-1960" "TMe-2004" "TMe-1906"
#> [37] "TMe-670"  "TMe-3142" "TMe-3463" "TMe-1272" "TMe-3462" "TMe-3264"
#> [43] "TMe-501"  "TMe-3219" "TMe-2383" "TMe-3351" "TMe-2131" "TMe-3475"
#> [49] "TMe-1224" "TMe-1218" "TMe-3394" "TMe-3132" "TMe-3553" "TMe-1466"
#> [55] "TMe-3398" "TMe-1140" "TMe-117"  "TMe-469"  "TMe-3623" "TMe-3479"
#> [61] "TMe-3416"
#> 
#> $II
#>  [1] "TMe-2414" "TMe-1732" "TMe-636"  "TMe-2862" "TMe-2952" "TMe-768" 
#>  [7] "TMe-1184" "TMe-3187" "TMe-3557" "TMe-930"  "TMe-1459" "TMe-2903"
#> [13] "TMe-86"   "TMe-3467" "TMe-2860" "TMe-2033" "TMe-1616" "TMe-2056"
#> [19] "TMe-2352" "TMe-171"  "TMe-2258" "TMe-2765" "TMe-3101" "TMe-1251"
#> [25] "TMe-1698" "TMe-2166" "TMe-3338" "TMe-3401" "TMe-3020" "TMe-3371"
#> [31] "TMe-431"  "TMe-3237" "TMe-509"  "TMe-1474" "TMe-1795" "TMe-3066"
#> [37] "TMe-3612" "TMe-2978" "TMe-3530" "TMe-681"  "TMe-3605"
#> 
#> $III
#>  [1] "TMe-3592" "TMe-1956" "TMe-3174" "TMe-3230" "TMe-3346" "TMe-3326"
#>  [7] "TMe-3721" "TMe-3274" "TMe-3556" "TMe-64"   "TMe-2502" "TMe-261" 
#> [13] "TMe-773"  "TMe-1200" "TMe-420"  "TMe-1993" "TMe-2751" "TMe-3383"
#> [19] "TMe-3029" "TMe-1787" "TMe-1965" "TMe-2394" "TMe-1863" "TMe-2811"
#> [25] "TMe-3445" "TMe-3731" "TMe-3407" "TMe-4"    "TMe-1176" "TMe-3118"
#> [31] "TMe-3575" "TMe-3069" "TMe-3397" "TMe-405"  "TMe-3376" "TMe-3100"
#> [37] "TMe-1910"
#> 
#> $IV
#>  [1] "TMe-34"   "TMe-801"  "TMe-3494" "TMe-352"  "TMe-3257" "TMe-519" 
#>  [7] "TMe-2775" "TMe-1267" "TMe-170"  "TMe-1987" "TMe-15"   "TMe-2332"
#> [13] "TMe-1297" "TMe-3511" "TMe-3779" "TMe-1700" "TMe-1151" "TMe-18"  
#> [19] "TMe-684"  "TMe-25"   "TMe-3089" "TMe-698"  "TMe-1579" "TMe-270" 
#> [25] "TMe-3660" "TMe-415"  "TMe-513"  "TMe-154"  "TMe-3760" "TMe-2039"
#> [31] "TMe-1020" "TMe-1010" "TMe-3440" "TMe-812"  "TMe-708"  "TMe-3017"
#> [37] "TMe-2924" "TMe-656"  "TMe-383"  "TMe-286"  "TMe-1479" "TMe-3435"
#> [43] "TMe-714"  "TMe-1903" "TMe-3639" "TMe-2390" "TMe-3558" "TMe-54"  
#> [49] "TMe-1078" "TMe-380"  "TMe-641"  "TMe-1330" "TMe-153"  "TMe-3276"
#> [55] "TMe-3757" "TMe-760"  "TMe-1553" "TMe-1278" "TMe-1129" "TMe-2764"
#> [61] "TMe-3278" "TMe-2210" "TMe-186"  "TMe-1577" "TMe-3591" "TMe-3225"
#> [67] "TMe-1397" "TMe-3541" "TMe-386"  "TMe-107"  "TMe-619"  "TMe-3054"
#> [73] "TMe-516"  "TMe-81"   "TMe-209"  "TMe-3084" "TMe-2971" "TMe-3340"
#> [79] "TMe-2981" "TMe-3161" "TMe-372" 
#> 
#> $V
#>  [1] "TMe-2018" "TMe-551"  "TMe-2032" "TMe-645"  "TMe-158"  "TMe-1299"
#>  [7] "TMe-323"  "TMe-3408" "TMe-1388" "TMe-390"  "TMe-1694" "TMe-547" 
#> [13] "TMe-245"  "TMe-929"  "TMe-472"  "TMe-1877" "TMe-2138" "TMe-1629"
#> [19] "TMe-2855" "TMe-1762" "TMe-2307" "TMe-627"  "TMe-1098" "TMe-1268"
#> [25] "TMe-2904" "TMe-1227" "TMe-2953" "TMe-194"  "TMe-1522" "TMe-892" 
#> [31] "TMe-723"  "TMe-2319" "TMe-1834" "TMe-385"  "TMe-1488" "TMe-2057"
#> [37] "TMe-2688" "TMe-795"  "TMe-771"  "TMe-2961" "TMe-2853" "TMe-370" 
#> [43] "TMe-603"  "TMe-828"  "TMe-424"  "TMe-901"  "TMe-687"  "TMe-1446"
#> [49] "TMe-1269" "TMe-475"  "TMe-1367" "TMe-826"  "TMe-1199" "TMe-943" 
#> [55] "TMe-287"  "TMe-870"  "TMe-1487" "TMe-755"  "TMe-924"  "TMe-326" 
#> [61] "TMe-825"  "TMe-2192" "TMe-969"  "TMe-1012" "TMe-1427" "TMe-976" 
#> [67] "TMe-439"  "TMe-432"  "TMe-1788" "TMe-1357" "TMe-418"  "TMe-884" 
#> [73] "TMe-823"  "TMe-193"  "TMe-1778" "TMe-1037" "TMe-1556" "TMe-412" 
#> [79] "TMe-3356" "TMe-2296"
#> 
#> $VI
#>  [1] "TMe-856"  "TMe-1816" "TMe-846"  "TMe-907"  "TMe-2067" "TMe-236" 
#>  [7] "TMe-620"  "TMe-2402" "TMe-1835" "TMe-3387" "TMe-893"  "TMe-836" 
#> [13] "TMe-1286" "TMe-1509" "TMe-2963" "TMe-213"  "TMe-1216" "TMe-2064"
#> [19] "TMe-2983" "TMe-662"  "TMe-876"  "TMe-514"  "TMe-728"  "TMe-693" 
#> [25] "TMe-1076" "TMe-1506" "TMe-2510" "TMe-985"  "TMe-3177" "TMe-752" 
#> [31] "TMe-3095" "TMe-631"  "TMe-531"  "TMe-1809" "TMe-1416" "TMe-854" 
#> [37] "TMe-1178"
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
