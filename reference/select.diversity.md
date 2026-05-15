# Selection of Entries from Clusters/Groups on the basis of Optimized Diversity

Select entries from cluster/groups in the entire collection which form a
subset with the highest trait diversity according to a either pooled or
mean diversity index estimate.

## Usage

``` r
select.diversity(
  data,
  names,
  group,
  alloc,
  qualitative,
  always.selected = NULL,
  div.index = c("richness", "shannon", "simpson", "mcintosh"),
  shannon.base = exp(1),
  div.fun = NULL,
  metric = c("mean", "pooled"),
  search = c("random", "greedy"),
  local.search = c("best.improvement", "first.improvement"),
  n.iter = 1000,
  max.iter = 30
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

- alloc:

  A named numeric vector specifying the number of entries to be
  selected. Names should correspond to the levels of the "`"group"`
  column, and values indicate the number of elements to be selected from
  each level.

- qualitative:

  Name of columns with the qualitative traits as a character vector.

- always.selected:

  Names of accessions to be always included in the core set as a
  character vector.

- div.index:

  The diversity index to be used to estimate within cluster/group
  diversity.

- shannon.base:

  The logarithm base to be used for estimation of Shannon diversity
  index. Default is `exp(1)`.

- div.fun:

  A function to estimate diversity index from a factor vector of
  qualitative trait data.

- metric:

  The metric to be computed from the diversity index. Either `"pooled"`
  or `"mean"`.

- search:

  Character string specifying the search strategy used to find the
  subset with the highest diversity score. Either `"random"` (default)
  or `"greedy"` (See **Details**).

- local.search:

  Character string specifying the local search strategy used in the
  1-opt improvement phase of the greedy search (`search = "greedy"`).
  Either `"best.improvement"` (default) or `"first.improvement"`.
  Ignored when `search = "random"`.

- n.iter:

  Integer specifying the number of random candidate subsets generated
  per group to optimze the diversity for random search
  (`search = "random"`).

- max.iter:

  The maximum number of 1-opt passes for greedy search
  (`search = "greedy"`).

## Value

A named list where each element contains the selected entry identifiers
for a cluster/group.

## Details

To identify subsets with highest diversity estimates, the following
strategies are available. These strategies are similar to the
"Maximization" or M strategy of Schoen and Brown (1993) .

### Random search / Monte Carlo Method

For each cluster/group, multiple candidate subsets are sampled randomly
and the subset with the highest trait diversity according to either
pooled or mean diversity index estimate is retained. The quality of the
solution improves with increasing `n.iter` but is not guaranteed to find
the global optimum (Anatoly Zhigljavsky and Antanas Zilinskas 2008) .

### Greedy search with 1-opt

This method builds a solution incrementally by adding the accession that
maximises the diversity score at each step, starting from the
`always.selected` accessions (or a single randomly drawn accession when
there are no accessions specified in `always.selected`) present in the
particular cluster/group (Nemhauser et al. 1978; Fisher et al. 1978;
Cormen et al. 2022) . The 'greedy' solution is then refined by a 1-opt
local search controlled by `local.search` and `max.iter` (Lin 1965) .
Greedy search is deterministic given a fixed `always.selected` set; when
there are no accessions specified in `always.selected` present in the
particular cluster/group results may vary across runs due to the random
initialisation.

`local.search = "best.improvement"` scans all possible single swaps in
each pass and applies the one yielding the greatest improvement before
restarting. his guarantees the steepest ascent at each pass but requires
evaluating all \\k \times (n - k)\\ swap pairs per pass, where \\k\\ is
the number of swappable accessions and \\n - k\\ is the size of the
candidate pool (Papadimitriou and Steiglitz 1998) .

`local.search = "first.improvement"` applies the first swap that
improves the score and immediately restarts the search. This typically
requires fewer score evaluations per pass and converges faster, but may
find a different local optimum than `"best.improvement"` (Papadimitriou
and Steiglitz 1998) .

Both strategies terminate when no improving swap exists (local optimum)
or when `max.iter` passes have been completed.

Entries listed as `always.selected` are mandatorily included in the
selection. Warnings are issued if requested allocation is smaller than
the number of always-selected entries in a cluster/group and/or when the
cluster/group does not contain enough remaining entries to fulfill the
allocation.

## References

Anatoly Zhigljavsky, Antanas Zilinskas (2008). *Stochastic Global
Optimization*, volume 9 of *Springer Optimization and Its Applications*.
Springer US, Boston, MA. ISBN 978-0-387-74022-5.  
  
Cormen TH, Leiserson CE, Rivest RL, Stein C (2022). *Introduction to
Algorithms*, 4 edition. MIT Press, Cambridge, MA, USA. ISBN
978-0-262-04630-5.  
  
Fisher ML, Nemhauser GL, Wolsey LA (1978). “An analysis of
approximations for maximizing submodular set functions-II.”
*Mathematical Programming Study*, **8**, 73–87.  
  
Lin S (1965). “Computer solutions of the traveling salesman problem.”
*Bell System Technical Journal*, **44**(10), 2245–2269.  
  
Nemhauser GL, Wolsey LA, Fisher ML (1978). “An analysis of
approximations for maximizing submodular set functions-I.” *Mathematical
Programming*, **14**(1), 265–294.  
  
Papadimitriou CH, Steiglitz K (1998). *Combinatorial optimization:
Algorithms and complexity*. Dover Publications, Mineola, N.Y. ISBN
978-0-486-40258-1.  
  
Schoen DJ, Brown AHD (1993). “Conservation of allelic richness in wild
crop relatives is aided by assessment of genetic markers.” *Proceedings
of the National Academy of Sciences*, **90**(22), 10623–10627.

## See also

[`select.random`](https://aravind-j.github.io/SampleCore/reference/select.random.md),
[`select.distance`](https://aravind-j.github.io/SampleCore/reference/select.distance.md)

## Examples

``` r

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare example data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(cluster)
library(ggplot2)

data(cassava_EC_gp)

data <- cbind(genotypes = rownames(cassava_EC_gp), cassava_EC_gp)

quant <- c("NMSR", "TTRN", "TFWSR", "TTRW", "TFWSS", "TTSW", "TTPW", "AVPW",
           "ARSR", "SRDM")
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

# Prepare inputs
counts <- c(I = 31, II = 31, III = 18, IV = 35, V = 40, VI = 17)

mand_accns <-
  c("TMe-34", "TMe-3423", "TMe-2018", "TMe-801", "TMe-551")

# Get distance matrix - Only for visualization

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Custom Diversity functions
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
# Random search
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mean richness
randomsel_mean_richness <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "richness",
                   metric = "mean", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_mean_richness
#> [[1]]
#>  [1] "TMe-3423" "TMe-2084" "TMe-3460" "TMe-2916" "TMe-2993" "TMe-3437"
#>  [7] "TMe-3698" "TMe-2004" "TMe-583"  "TMe-1906" "TMe-300"  "TMe-1915"
#> [13] "TMe-3477" "TMe-2779" "TMe-3534" "TMe-1451" "TMe-3475" "TMe-1716"
#> [19] "TMe-1922" "TMe-2069" "TMe-1218" "TMe-642"  "TMe-3389" "TMe-2453"
#> [25] "TMe-2810" "TMe-1224" "TMe-3462" "TMe-407"  "TMe-3694" "TMe-1360"
#> [31] "TMe-1425"
#> 
#> [[2]]
#>  [1] "TMe-74"   "TMe-160"  "TMe-3455" "TMe-1739" "TMe-3146" "TMe-2412"
#>  [7] "TMe-2970" "TMe-2998" "TMe-171"  "TMe-1831" "TMe-3101" "TMe-2851"
#> [13] "TMe-1184" "TMe-2705" "TMe-3547" "TMe-3690" "TMe-1668" "TMe-3564"
#> [19] "TMe-601"  "TMe-3540" "TMe-3286" "TMe-369"  "TMe-3021" "TMe-3140"
#> [25] "TMe-409"  "TMe-3466" "TMe-2952" "TMe-3308" "TMe-664"  "TMe-2860"
#> [31] "TMe-2891"
#> 
#> [[3]]
#>  [1] "TMe-946"  "TMe-261"  "TMe-3638" "TMe-3575" "TMe-3643" "TMe-234" 
#>  [7] "TMe-64"   "TMe-3663" "TMe-3148" "TMe-3362" "TMe-1965" "TMe-3008"
#> [13] "TMe-1897" "TMe-2356" "TMe-434"  "TMe-3335" "TMe-13"   "TMe-2161"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-3253" "TMe-352"  "TMe-1434" "TMe-241" 
#>  [7] "TMe-2924" "TMe-226"  "TMe-2020" "TMe-1369" "TMe-2839" "TMe-43"  
#> [13] "TMe-2340" "TMe-1988" "TMe-1078" "TMe-735"  "TMe-3327" "TMe-3204"
#> [19] "TMe-386"  "TMe-3109" "TMe-3238" "TMe-1027" "TMe-3214" "TMe-403" 
#> [25] "TMe-107"  "TMe-1297" "TMe-480"  "TMe-2809" "TMe-2059" "TMe-1150"
#> [31] "TMe-2960" "TMe-594"  "TMe-1419" "TMe-1123" "TMe-761" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2853" "TMe-1879" "TMe-1308" "TMe-3408"
#>  [7] "TMe-1375" "TMe-3185" "TMe-1609" "TMe-2790" "TMe-868"  "TMe-189" 
#> [13] "TMe-673"  "TMe-378"  "TMe-472"  "TMe-423"  "TMe-2753" "TMe-1418"
#> [19] "TMe-764"  "TMe-1268" "TMe-2057" "TMe-2612" "TMe-863"  "TMe-889" 
#> [25] "TMe-886"  "TMe-3034" "TMe-162"  "TMe-688"  "TMe-629"  "TMe-382" 
#> [31] "TMe-2905" "TMe-412"  "TMe-474"  "TMe-826"  "TMe-2192" "TMe-2060"
#> [37] "TMe-2590" "TMe-565"  "TMe-2542" "TMe-1401"
#> 
#> [[6]]
#>  [1] "TMe-3549" "TMe-315"  "TMe-185"  "TMe-1383" "TMe-470"  "TMe-1483"
#>  [7] "TMe-1264" "TMe-465"  "TMe-1875" "TMe-690"  "TMe-936"  "TMe-659" 
#> [13] "TMe-1413" "TMe-2551" "TMe-1539" "TMe-696"  "TMe-1558"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search", subtitle = "Mean richness")


# Pooled richness
randomsel_sum_richness <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "richness",
                   metric = "pooled", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_sum_richness
#> [[1]]
#>  [1] "TMe-3423" "TMe-3515" "TMe-716"  "TMe-3002" "TMe-3003" "TMe-2453"
#>  [7] "TMe-3625" "TMe-2031" "TMe-3030" "TMe-1717" "TMe-1272" "TMe-299" 
#> [13] "TMe-3169" "TMe-1224" "TMe-3145" "TMe-1360" "TMe-990"  "TMe-3553"
#> [19] "TMe-757"  "TMe-3074" "TMe-1425" "TMe-3493" "TMe-3223" "TMe-3441"
#> [25] "TMe-2984" "TMe-1786" "TMe-2985" "TMe-3149" "TMe-469"  "TMe-2965"
#> [31] "TMe-2773"
#> 
#> [[2]]
#>  [1] "TMe-2048" "TMe-3141" "TMe-3805" "TMe-47"   "TMe-2998" "TMe-664" 
#>  [7] "TMe-3210" "TMe-2056" "TMe-1833" "TMe-1831" "TMe-3338" "TMe-1242"
#> [13] "TMe-3466" "TMe-2412" "TMe-1271" "TMe-171"  "TMe-74"   "TMe-1795"
#> [19] "TMe-3150" "TMe-3605" "TMe-3800" "TMe-2915" "TMe-3101" "TMe-59"  
#> [25] "TMe-3286" "TMe-3557" "TMe-509"  "TMe-3766" "TMe-3371" "TMe-2166"
#> [31] "TMe-2258"
#> 
#> [[3]]
#>  [1] "TMe-3715" "TMe-3336" "TMe-3274" "TMe-2756" "TMe-1398" "TMe-3005"
#>  [7] "TMe-187"  "TMe-1288" "TMe-773"  "TMe-3287" "TMe-946"  "TMe-3207"
#> [13] "TMe-617"  "TMe-2502" "TMe-2977" "TMe-35"   "TMe-3234" "TMe-3556"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-3409" "TMe-3545" "TMe-3290" "TMe-3212"
#>  [7] "TMe-2059" "TMe-2841" "TMe-394"  "TMe-1297" "TMe-3068" "TMe-1485"
#> [13] "TMe-1027" "TMe-3494" "TMe-241"  "TMe-2337" "TMe-3511" "TMe-3390"
#> [19] "TMe-2378" "TMe-1348" "TMe-1579" "TMe-3442" "TMe-2922" "TMe-3089"
#> [25] "TMe-184"  "TMe-1139" "TMe-1267" "TMe-2971" "TMe-1368" "TMe-225" 
#> [31] "TMe-1525" "TMe-3707" "TMe-733"  "TMe-352"  "TMe-594" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1411" "TMe-1760" "TMe-1004" "TMe-476" 
#>  [7] "TMe-347"  "TMe-399"  "TMe-212"  "TMe-412"  "TMe-2355" "TMe-457" 
#> [13] "TMe-1293" "TMe-712"  "TMe-333"  "TMe-142"  "TMe-479"  "TMe-1957"
#> [19] "TMe-1257" "TMe-3185" "TMe-307"  "TMe-600"  "TMe-655"  "TMe-1268"
#> [25] "TMe-1934" "TMe-211"  "TMe-2769" "TMe-1227" "TMe-588"  "TMe-363" 
#> [31] "TMe-1547" "TMe-1446" "TMe-1388" "TMe-948"  "TMe-1312" "TMe-2542"
#> [37] "TMe-764"  "TMe-224"  "TMe-3628" "TMe-742" 
#> 
#> [[6]]
#>  [1] "TMe-726"  "TMe-1256" "TMe-690"  "TMe-2195" "TMe-222"  "TMe-1174"
#>  [7] "TMe-1608" "TMe-963"  "TMe-1472" "TMe-1775" "TMe-662"  "TMe-1428"
#> [13] "TMe-696"  "TMe-315"  "TMe-693"  "TMe-1847" "TMe-268" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search", subtitle = "Pooled richness")


# Mean Shannon-Weaver diversity index
randomsel_mean_shannon <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "shannon",
                   metric = "mean", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_mean_shannon
#> [[1]]
#>  [1] "TMe-3423" "TMe-3534" "TMe-1935" "TMe-3132" "TMe-1247" "TMe-3236"
#>  [7] "TMe-3475" "TMe-606"  "TMe-1272" "TMe-2993" "TMe-3110" "TMe-2010"
#> [13] "TMe-1931" "TMe-2152" "TMe-3490" "TMe-3164" "TMe-1964" "TMe-486" 
#> [19] "TMe-3111" "TMe-3457" "TMe-3389" "TMe-3337" "TMe-264"  "TMe-3292"
#> [25] "TMe-1914" "TMe-2066" "TMe-756"  "TMe-2966" "TMe-3359" "TMe-41"  
#> [31] "TMe-3462"
#> 
#> [[2]]
#>  [1] "TMe-2021" "TMe-2056" "TMe-3530" "TMe-1251" "TMe-960"  "TMe-2935"
#>  [7] "TMe-509"  "TMe-160"  "TMe-2952" "TMe-3239" "TMe-2995" "TMe-3286"
#> [13] "TMe-842"  "TMe-2862" "TMe-3140" "TMe-3237" "TMe-3466" "TMe-1005"
#> [19] "TMe-2978" "TMe-2166" "TMe-3200" "TMe-3210" "TMe-1732" "TMe-289" 
#> [25] "TMe-1172" "TMe-1504" "TMe-2860" "TMe-3771" "TMe-339"  "TMe-636" 
#> [31] "TMe-3338"
#> 
#> [[3]]
#>  [1] "TMe-1262" "TMe-3397" "TMe-3143" "TMe-2405" "TMe-3346" "TMe-1176"
#>  [7] "TMe-1897" "TMe-1819" "TMe-2977" "TMe-3551" "TMe-785"  "TMe-3565"
#> [13] "TMe-3326" "TMe-3804" "TMe-2756" "TMe-3230" "TMe-2270" "TMe-3287"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-1456" "TMe-150"  "TMe-3276" "TMe-3707"
#>  [7] "TMe-1010" "TMe-38"   "TMe-225"  "TMe-1525" "TMe-3278" "TMe-27"  
#> [13] "TMe-3545" "TMe-513"  "TMe-1336" "TMe-1348" "TMe-1828" "TMe-82"  
#> [19] "TMe-656"  "TMe-3390" "TMe-3538" "TMe-467"  "TMe-3269" "TMe-3004"
#> [25] "TMe-280"  "TMe-2025" "TMe-2240" "TMe-3168" "TMe-3378" "TMe-3144"
#> [31] "TMe-3054" "TMe-3000" "TMe-1162" "TMe-3212" "TMe-3760"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-797"  "TMe-1158" "TMe-823"  "TMe-1879"
#>  [7] "TMe-2589" "TMe-1521" "TMe-2413" "TMe-932"  "TMe-360"  "TMe-2959"
#> [13] "TMe-1411" "TMe-584"  "TMe-627"  "TMe-858"  "TMe-277"  "TMe-48"  
#> [19] "TMe-782"  "TMe-3411" "TMe-476"  "TMe-1834" "TMe-137"  "TMe-1488"
#> [25] "TMe-1382" "TMe-1183" "TMe-2055" "TMe-731"  "TMe-2820" "TMe-1622"
#> [31] "TMe-954"  "TMe-1862" "TMe-1265" "TMe-3408" "TMe-3185" "TMe-1455"
#> [37] "TMe-651"  "TMe-1523" "TMe-3356" "TMe-99"  
#> 
#> [[6]]
#>  [1] "TMe-1531" "TMe-1483" "TMe-1832" "TMe-315"  "TMe-1302" "TMe-1124"
#>  [7] "TMe-1518" "TMe-1566" "TMe-856"  "TMe-883"  "TMe-1239" "TMe-1238"
#> [13] "TMe-269"  "TMe-2064" "TMe-1445" "TMe-2543" "TMe-1509"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Mean Shannon-Weaver diversity index")


# Pooled Shannon-Weaver diversity index
randomsel_sum_shannon <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "shannon",
                   metric = "pooled", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_sum_shannon
#> [[1]]
#>  [1] "TMe-3423" "TMe-2909" "TMe-1964" "TMe-2066" "TMe-3026" "TMe-2993"
#>  [7] "TMe-2069" "TMe-1425" "TMe-1930" "TMe-841"  "TMe-3724" "TMe-1247"
#> [13] "TMe-3429" "TMe-3396" "TMe-3325" "TMe-3229" "TMe-717"  "TMe-3252"
#> [19] "TMe-606"  "TMe-2927" "TMe-3264" "TMe-2779" "TMe-1886" "TMe-2975"
#> [25] "TMe-3481" "TMe-290"  "TMe-2453" "TMe-3490" "TMe-579"  "TMe-2785"
#> [31] "TMe-1869"
#> 
#> [[2]]
#>  [1] "TMe-902"  "TMe-3338" "TMe-1444" "TMe-1005" "TMe-1754" "TMe-1385"
#>  [7] "TMe-1919" "TMe-3210" "TMe-2054" "TMe-2851" "TMe-1831" "TMe-2757"
#> [13] "TMe-2352" "TMe-2211" "TMe-2997" "TMe-2428" "TMe-86"   "TMe-1732"
#> [19] "TMe-1323" "TMe-2568" "TMe-2998" "TMe-1619" "TMe-1505" "TMe-2866"
#> [25] "TMe-2935" "TMe-3766" "TMe-2715" "TMe-369"  "TMe-3200" "TMe-2021"
#> [31] "TMe-2033"
#> 
#> [[3]]
#>  [1] "TMe-785"  "TMe-2756" "TMe-1443" "TMe-3458" "TMe-1176" "TMe-3556"
#>  [7] "TMe-3565" "TMe-3383" "TMe-174"  "TMe-1593" "TMe-3100" "TMe-3133"
#> [13] "TMe-187"  "TMe-1804" "TMe-234"  "TMe-1738" "TMe-2169" "TMe-635" 
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-3054" "TMe-63"   "TMe-103" 
#>  [7] "TMe-967"  "TMe-210"  "TMe-372"  "TMe-3758" "TMe-1406" "TMe-66"  
#> [13] "TMe-1479" "TMe-891"  "TMe-734"  "TMe-1123" "TMe-24"   "TMe-3025"
#> [19] "TMe-241"  "TMe-3040" "TMe-2020" "TMe-897"  "TMe-18"   "TMe-1376"
#> [25] "TMe-2988" "TMe-3167" "TMe-467"  "TMe-2924" "TMe-3290" "TMe-824" 
#> [31] "TMe-650"  "TMe-7"    "TMe-2979" "TMe-3269" "TMe-235" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-432"  "TMe-804"  "TMe-536"  "TMe-217" 
#>  [7] "TMe-2750" "TMe-2612" "TMe-723"  "TMe-2855" "TMe-2003" "TMe-194" 
#> [13] "TMe-2060" "TMe-798"  "TMe-2590" "TMe-2435" "TMe-600"  "TMe-764" 
#> [19] "TMe-1812" "TMe-1311" "TMe-777"  "TMe-1488" "TMe-816"  "TMe-1523"
#> [25] "TMe-1500" "TMe-1099" "TMe-2016" "TMe-822"  "TMe-373"  "TMe-819" 
#> [31] "TMe-385"  "TMe-1227" "TMe-162"  "TMe-1188" "TMe-782"  "TMe-669" 
#> [37] "TMe-2413" "TMe-2863" "TMe-748"  "TMe-343" 
#> 
#> [[6]]
#>  [1] "TMe-1261" "TMe-1033" "TMe-3095" "TMe-1518" "TMe-1178" "TMe-845" 
#>  [7] "TMe-1548" "TMe-1182" "TMe-1472" "TMe-693"  "TMe-3177" "TMe-185" 
#> [13] "TMe-2510" "TMe-983"  "TMe-626"  "TMe-310"  "TMe-1678"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Pooled Shannon-Weaver diversity index")


# Mean Gini-Simpson diversity index
randomsel_mean_simpson <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "simpson",
                   metric = "mean", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_mean_simpson
#> [[1]]
#>  [1] "TMe-3423" "TMe-2976" "TMe-3169" "TMe-3111" "TMe-2912" "TMe-3515"
#>  [7] "TMe-3235" "TMe-839"  "TMe-1581" "TMe-3719" "TMe-3518" "TMe-1734"
#> [13] "TMe-3726" "TMe-248"  "TMe-3319" "TMe-2069" "TMe-3433" "TMe-815" 
#> [19] "TMe-2993" "TMe-1960" "TMe-1958" "TMe-3208" "TMe-3252" "TMe-3236"
#> [25] "TMe-3694" "TMe-1272" "TMe-3479" "TMe-2984" "TMe-2513" "TMe-3292"
#> [31] "TMe-3323"
#> 
#> [[2]]
#>  [1] "TMe-3141" "TMe-3209" "TMe-648"  "TMe-3308" "TMe-3467" "TMe-3009"
#>  [7] "TMe-3401" "TMe-176"  "TMe-3771" "TMe-2414" "TMe-2257" "TMe-636" 
#> [13] "TMe-2128" "TMe-3690" "TMe-47"   "TMe-1323" "TMe-1833" "TMe-2123"
#> [19] "TMe-3605" "TMe-1505" "TMe-3533" "TMe-1831" "TMe-902"  "TMe-2851"
#> [25] "TMe-3188" "TMe-2860" "TMe-539"  "TMe-1107" "TMe-171"  "TMe-3237"
#> [31] "TMe-1005"
#> 
#> [[3]]
#>  [1] "TMe-3700" "TMe-2362" "TMe-3663" "TMe-187"  "TMe-1790" "TMe-2285"
#>  [7] "TMe-381"  "TMe-2154" "TMe-3128" "TMe-3596" "TMe-261"  "TMe-1804"
#> [13] "TMe-3299" "TMe-14"   "TMe-1897" "TMe-1838" "TMe-617"  "TMe-1787"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-766"  "TMe-684"  "TMe-3511" "TMe-2567"
#>  [7] "TMe-550"  "TMe-1016" "TMe-3576" "TMe-1397" "TMe-3757" "TMe-2776"
#> [13] "TMe-2025" "TMe-107"  "TMe-3760" "TMe-3044" "TMe-3602" "TMe-1434"
#> [19] "TMe-2337" "TMe-2210" "TMe-24"   "TMe-3025" "TMe-2958" "TMe-317" 
#> [25] "TMe-794"  "TMe-2979" "TMe-3459" "TMe-3167" "TMe-919"  "TMe-1402"
#> [31] "TMe-352"  "TMe-3072" "TMe-1231" "TMe-388"  "TMe-3537"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1802" "TMe-359"  "TMe-1534" "TMe-755" 
#>  [7] "TMe-1305" "TMe-284"  "TMe-1411" "TMe-424"  "TMe-1873" "TMe-1257"
#> [13] "TMe-334"  "TMe-673"  "TMe-2907" "TMe-2041" "TMe-2413" "TMe-651" 
#> [19] "TMe-1248" "TMe-1871" "TMe-448"  "TMe-700"  "TMe-2271" "TMe-627" 
#> [25] "TMe-323"  "TMe-217"  "TMe-1373" "TMe-584"  "TMe-803"  "TMe-2121"
#> [31] "TMe-194"  "TMe-326"  "TMe-2057" "TMe-224"  "TMe-679"  "TMe-399" 
#> [37] "TMe-423"  "TMe-307"  "TMe-2853" "TMe-3411"
#> 
#> [[6]]
#>  [1] "TMe-1182" "TMe-2534" "TMe-1992" "TMe-2572" "TMe-1661" "TMe-1902"
#>  [7] "TMe-1124" "TMe-213"  "TMe-531"  "TMe-1445" "TMe-2543" "TMe-1383"
#> [13] "TMe-2026" "TMe-470"  "TMe-315"  "TMe-2035" "TMe-661" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Mean Gini-Simpson diversity index")


# Pooled Gini-Simpson diversity index
randomsel_sum_simpson <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "simpson",
                   metric = "pooled", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_sum_simpson
#> [[1]]
#>  [1] "TMe-3423" "TMe-3236" "TMe-677"  "TMe-2967" "TMe-3415" "TMe-3437"
#>  [7] "TMe-3266" "TMe-2996" "TMe-3112" "TMe-3229" "TMe-727"  "TMe-290" 
#> [13] "TMe-3462" "TMe-756"  "TMe-2031" "TMe-1341" "TMe-2004" "TMe-1451"
#> [19] "TMe-3345" "TMe-41"   "TMe-2206" "TMe-1190" "TMe-2943" "TMe-1347"
#> [25] "TMe-3296" "TMe-3705" "TMe-2103" "TMe-3485" "TMe-3429" "TMe-2949"
#> [31] "TMe-3353"
#> 
#> [[2]]
#>  [1] "TMe-1907" "TMe-636"  "TMe-842"  "TMe-648"  "TMe-2033" "TMe-2166"
#>  [7] "TMe-1860" "TMe-39"   "TMe-3338" "TMe-53"   "TMe-1969" "TMe-3766"
#> [13] "TMe-2851" "TMe-339"  "TMe-3209" "TMe-3642" "TMe-3188" "TMe-2891"
#> [19] "TMe-2970" "TMe-2412" "TMe-2950" "TMe-47"   "TMe-1668" "TMe-3020"
#> [25] "TMe-1833" "TMe-3466" "TMe-2054" "TMe-2058" "TMe-2268" "TMe-1385"
#> [31] "TMe-2329"
#> 
#> [[3]]
#>  [1] "TMe-381"  "TMe-3174" "TMe-2977" "TMe-1443" "TMe-2926" "TMe-3362"
#>  [7] "TMe-1792" "TMe-3631" "TMe-3397" "TMe-2154" "TMe-2968" "TMe-3663"
#> [13] "TMe-3383" "TMe-1897" "TMe-3572" "TMe-1868" "TMe-3046" "TMe-3551"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-2410" "TMe-1335" "TMe-1470" "TMe-3541"
#>  [7] "TMe-2980" "TMe-2947" "TMe-3055" "TMe-684"  "TMe-1267" "TMe-3591"
#> [13] "TMe-3581" "TMe-30"   "TMe-1020" "TMe-641"  "TMe-226"  "TMe-1139"
#> [19] "TMe-280"  "TMe-806"  "TMe-3168" "TMe-1348" "TMe-3273" "TMe-3044"
#> [25] "TMe-2924" "TMe-241"  "TMe-3054" "TMe-3204" "TMe-266"  "TMe-3602"
#> [31] "TMe-65"   "TMe-3615" "TMe-3242" "TMe-2788" "TMe-1278"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1192" "TMe-1012" "TMe-355"  "TMe-1399"
#>  [7] "TMe-1257" "TMe-574"  "TMe-347"  "TMe-2413" "TMe-391"  "TMe-1375"
#> [13] "TMe-1440" "TMe-889"  "TMe-1488" "TMe-1427" "TMe-2339" "TMe-343" 
#> [19] "TMe-1367" "TMe-1157" "TMe-1220" "TMe-193"  "TMe-1160" "TMe-1924"
#> [25] "TMe-224"  "TMe-2425" "TMe-954"  "TMe-800"  "TMe-2121" "TMe-2435"
#> [31] "TMe-723"  "TMe-1215" "TMe-822"  "TMe-1299" "TMe-1963" "TMe-1391"
#> [37] "TMe-423"  "TMe-1901" "TMe-1722" "TMe-2904"
#> 
#> [[6]]
#>  [1] "TMe-1503" "TMe-580"  "TMe-1661" "TMe-1035" "TMe-2510" "TMe-1518"
#>  [7] "TMe-2534" "TMe-3387" "TMe-1416" "TMe-1124" "TMe-514"  "TMe-1362"
#> [13] "TMe-838"  "TMe-1033" "TMe-781"  "TMe-1875" "TMe-2572"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Pooled Gini-Simpson diversity index")


# Mean McIntosh diversity index
randomsel_mean_mcintosh <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "mcintosh",
                   metric = "pooled", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_mean_mcintosh
#> [[1]]
#>  [1] "TMe-3423" "TMe-3465" "TMe-3416" "TMe-3223" "TMe-1922" "TMe-2069"
#>  [7] "TMe-815"  "TMe-1564" "TMe-1190" "TMe-670"  "TMe-3115" "TMe-888" 
#> [13] "TMe-3488" "TMe-2253" "TMe-3394" "TMe-727"  "TMe-3521" "TMe-132" 
#> [19] "TMe-2936" "TMe-778"  "TMe-1906" "TMe-898"  "TMe-1468" "TMe-1568"
#> [25] "TMe-3292" "TMe-1589" "TMe-33"   "TMe-3163" "TMe-3418" "TMe-3633"
#> [31] "TMe-2462"
#> 
#> [[2]]
#>  [1] "TMe-636"  "TMe-3599" "TMe-3498" "TMe-2166" "TMe-1051" "TMe-1864"
#>  [7] "TMe-2998" "TMe-3188" "TMe-3766" "TMe-774"  "TMe-3258" "TMe-1860"
#> [13] "TMe-2085" "TMe-3150" "TMe-3009" "TMe-1969" "TMe-1444" "TMe-3800"
#> [19] "TMe-1241" "TMe-2"    "TMe-1619" "TMe-3371" "TMe-1107" "TMe-1668"
#> [25] "TMe-2797" "TMe-2412" "TMe-2123" "TMe-3210" "TMe-2891" "TMe-3021"
#> [31] "TMe-478" 
#> 
#> [[3]]
#>  [1] "TMe-1262" "TMe-405"  "TMe-3598" "TMe-3143" "TMe-3088" "TMe-3382"
#>  [7] "TMe-2811" "TMe-3458" "TMe-3804" "TMe-3556" "TMe-3274" "TMe-3750"
#> [13] "TMe-1398" "TMe-3569" "TMe-1147" "TMe-946"  "TMe-2347" "TMe-3005"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-3040" "TMe-3573" "TMe-619"  "TMe-761" 
#>  [7] "TMe-2947" "TMe-3422" "TMe-967"  "TMe-1903" "TMe-43"   "TMe-885" 
#> [13] "TMe-1406" "TMe-2841" "TMe-3019" "TMe-3758" "TMe-1335" "TMe-1434"
#> [19] "TMe-2025" "TMe-1364" "TMe-318"  "TMe-1948" "TMe-3276" "TMe-3537"
#> [25] "TMe-170"  "TMe-421"  "TMe-513"  "TMe-2410" "TMe-806"  "TMe-1155"
#> [31] "TMe-3435" "TMe-1330" "TMe-2890" "TMe-3242" "TMe-280" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2192" "TMe-1285" "TMe-1345" "TMe-2904"
#>  [7] "TMe-1026" "TMe-1778" "TMe-2863" "TMe-2612" "TMe-2319" "TMe-1183"
#> [13] "TMe-256"  "TMe-1332" "TMe-181"  "TMe-436"  "TMe-892"  "TMe-769" 
#> [19] "TMe-742"  "TMe-588"  "TMe-1367" "TMe-1420" "TMe-1722" "TMe-245" 
#> [25] "TMe-1101" "TMe-464"  "TMe-3659" "TMe-344"  "TMe-3311" "TMe-1307"
#> [31] "TMe-536"  "TMe-1901" "TMe-1446" "TMe-972"  "TMe-2050" "TMe-1159"
#> [37] "TMe-3411" "TMe-1501" "TMe-1963" "TMe-3356"
#> 
#> [[6]]
#>  [1] "TMe-1945" "TMe-838"  "TMe-1428" "TMe-728"  "TMe-2543" "TMe-2534"
#>  [7] "TMe-268"  "TMe-3095" "TMe-856"  "TMe-1035" "TMe-442"  "TMe-1509"
#> [13] "TMe-465"  "TMe-1796" "TMe-1554" "TMe-1416" "TMe-1769"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Mean McIntosh diversity index")


# Pooled McIntosh diversity index
randomsel_sum_mcintosh <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "mcintosh",
                   metric = "pooled", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_sum_mcintosh
#> [[1]]
#>  [1] "TMe-3423" "TMe-1425" "TMe-2993" "TMe-3726" "TMe-1190" "TMe-756" 
#>  [7] "TMe-3337" "TMe-41"   "TMe-841"  "TMe-642"  "TMe-3719" "TMe-3625"
#> [13] "TMe-1621" "TMe-3462" "TMe-3026" "TMe-1869" "TMe-2937" "TMe-1851"
#> [19] "TMe-3437" "TMe-3480" "TMe-2913" "TMe-2191" "TMe-2237" "TMe-882" 
#> [25] "TMe-3264" "TMe-2955" "TMe-888"  "TMe-3235" "TMe-1940" "TMe-3325"
#> [31] "TMe-2785"
#> 
#> [[2]]
#>  [1] "TMe-2268" "TMe-2056" "TMe-3203" "TMe-1668" "TMe-2077" "TMe-2127"
#>  [7] "TMe-3564" "TMe-2866" "TMe-60"   "TMe-3771" "TMe-2128" "TMe-3401"
#> [13] "TMe-3605" "TMe-3141" "TMe-1616" "TMe-3406" "TMe-960"  "TMe-3366"
#> [19] "TMe-3021" "TMe-2903" "TMe-1409" "TMe-3188" "TMe-674"  "TMe-2715"
#> [25] "TMe-3533" "TMe-404"  "TMe-3766" "TMe-2891" "TMe-3222" "TMe-1505"
#> [31] "TMe-2329"
#> 
#> [[3]]
#>  [1] "TMe-3700" "TMe-64"   "TMe-1910" "TMe-1838" "TMe-3287" "TMe-267" 
#>  [7] "TMe-3007" "TMe-2167" "TMe-3048" "TMe-174"  "TMe-123"  "TMe-1200"
#> [13] "TMe-853"  "TMe-3005" "TMe-3317" "TMe-1868" "TMe-3376" "TMe-1262"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-3204" "TMe-2399" "TMe-3257" "TMe-3273"
#>  [7] "TMe-15"   "TMe-153"  "TMe-2059" "TMe-3054" "TMe-2971" "TMe-1867"
#> [13] "TMe-1959" "TMe-3044" "TMe-388"  "TMe-967"  "TMe-2807" "TMe-3442"
#> [19] "TMe-594"  "TMe-1828" "TMe-3212" "TMe-1267" "TMe-550"  "TMe-2043"
#> [25] "TMe-383"  "TMe-1903" "TMe-3297" "TMe-3225" "TMe-1010" "TMe-3760"
#> [31] "TMe-2980" "TMe-3243" "TMe-241"  "TMe-103"  "TMe-3542"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1293" "TMe-2192" "TMe-1871" "TMe-782" 
#>  [7] "TMe-2835" "TMe-2769" "TMe-1802" "TMe-1359" "TMe-828"  "TMe-612" 
#> [13] "TMe-1453" "TMe-1357" "TMe-417"  "TMe-376"  "TMe-536"  "TMe-1381"
#> [19] "TMe-275"  "TMe-1762" "TMe-2753" "TMe-294"  "TMe-391"  "TMe-2761"
#> [25] "TMe-2589" "TMe-1482" "TMe-977"  "TMe-385"  "TMe-2959" "TMe-2307"
#> [31] "TMe-355"  "TMe-618"  "TMe-786"  "TMe-402"  "TMe-798"  "TMe-3332"
#> [37] "TMe-1934" "TMe-1500" "TMe-1332" "TMe-1458"
#> 
#> [[6]]
#>  [1] "TMe-1428" "TMe-2534" "TMe-580"  "TMe-361"  "TMe-1174" "TMe-321" 
#>  [7] "TMe-659"  "TMe-2963" "TMe-1486" "TMe-1007" "TMe-2398" "TMe-548" 
#> [13] "TMe-1392" "TMe-983"  "TMe-1518" "TMe-514"  "TMe-138" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Pooled McIntosh diversity index")


# Mean Brillouin diversity index
randomsel_mean_brillouin <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_brillouin,
                   metric = "mean", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_mean_brillouin
#> [[1]]
#>  [1] "TMe-3423" "TMe-2945" "TMe-2934" "TMe-2984" "TMe-3548" "TMe-2010"
#>  [7] "TMe-717"  "TMe-3496" "TMe-1973" "TMe-937"  "TMe-3625" "TMe-2027"
#> [13] "TMe-1360" "TMe-1851" "TMe-2513" "TMe-1224" "TMe-3132" "TMe-3051"
#> [19] "TMe-1140" "TMe-3396" "TMe-1915" "TMe-3488" "TMe-839"  "TMe-583" 
#> [25] "TMe-3478" "TMe-1156" "TMe-1600" "TMe-2069" "TMe-2985" "TMe-2913"
#> [31] "TMe-3223"
#> 
#> [[2]]
#>  [1] "TMe-2866" "TMe-3771" "TMe-528"  "TMe-404"  "TMe-1484" "TMe-3605"
#>  [7] "TMe-539"  "TMe-1323" "TMe-2995" "TMe-3766" "TMe-2860" "TMe-3237"
#> [13] "TMe-2200" "TMe-2077" "TMe-1831" "TMe-3612" "TMe-60"   "TMe-1619"
#> [19] "TMe-477"  "TMe-160"  "TMe-2"    "TMe-2058" "TMe-1005" "TMe-2258"
#> [25] "TMe-400"  "TMe-3443" "TMe-2952" "TMe-2903" "TMe-2054" "TMe-455" 
#> [31] "TMe-3599"
#> 
#> [[3]]
#>  [1] "TMe-3596" "TMe-187"  "TMe-3565" "TMe-3008" "TMe-3085" "TMe-3277"
#>  [7] "TMe-995"  "TMe-3556" "TMe-3048" "TMe-1897" "TMe-2926" "TMe-1939"
#> [13] "TMe-3575" "TMe-3274" "TMe-3551" "TMe-6"    "TMe-1207" "TMe-3028"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-737"  "TMe-180"  "TMe-3573" "TMe-3242"
#>  [7] "TMe-2971" "TMe-2924" "TMe-3255" "TMe-266"  "TMe-65"   "TMe-3055"
#> [13] "TMe-2928" "TMe-396"  "TMe-3538" "TMe-2809" "TMe-25"   "TMe-3006"
#> [19] "TMe-63"   "TMe-2947" "TMe-3256" "TMe-3535" "TMe-733"  "TMe-887" 
#> [25] "TMe-1364" "TMe-734"  "TMe-3269" "TMe-1891" "TMe-1027" "TMe-708" 
#> [31] "TMe-3576" "TMe-103"  "TMe-550"  "TMe-1330" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1273" "TMe-1963" "TMe-307"  "TMe-334" 
#>  [7] "TMe-1099" "TMe-612"  "TMe-456"  "TMe-2016" "TMe-98"   "TMe-344" 
#> [13] "TMe-312"  "TMe-943"  "TMe-1290" "TMe-1834" "TMe-275"  "TMe-1873"
#> [19] "TMe-1098" "TMe-1934" "TMe-3311" "TMe-2790" "TMe-527"  "TMe-2041"
#> [25] "TMe-355"  "TMe-1192" "TMe-336"  "TMe-647"  "TMe-1572" "TMe-2953"
#> [31] "TMe-2820" "TMe-2003" "TMe-723"  "TMe-168"  "TMe-1455" "TMe-1877"
#> [37] "TMe-972"  "TMe-376"  "TMe-2213" "TMe-3332"
#> 
#> [[6]]
#>  [1] "TMe-2983" "TMe-1403" "TMe-1174" "TMe-1566" "TMe-1902" "TMe-315" 
#>  [7] "TMe-1721" "TMe-659"  "TMe-836"  "TMe-2543" "TMe-1025" "TMe-728" 
#> [13] "TMe-726"  "TMe-1832" "TMe-2064" "TMe-531"  "TMe-2067"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Mean Brillouin diversity index")


# Pooled Brillouin diversity index
randomsel_sum_brillouin <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_brillouin,
                   metric = "pooled", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_sum_brillouin
#> [[1]]
#>  [1] "TMe-3423" "TMe-2993" "TMe-686"  "TMe-3396" "TMe-3481" "TMe-1096"
#>  [7] "TMe-1621" "TMe-2975" "TMe-3051" "TMe-41"   "TMe-2084" "TMe-3102"
#> [13] "TMe-2939" "TMe-3211" "TMe-3351" "TMe-3521" "TMe-407"  "TMe-2917"
#> [19] "TMe-3437" "TMe-1468" "TMe-3027" "TMe-3641" "TMe-1915" "TMe-3031"
#> [25] "TMe-3314" "TMe-3392" "TMe-1960" "TMe-756"  "TMe-2372" "TMe-2910"
#> [31] "TMe-1247"
#> 
#> [[2]]
#>  [1] "TMe-3533" "TMe-1137" "TMe-1907" "TMe-3495" "TMe-3114" "TMe-648" 
#>  [7] "TMe-2048" "TMe-60"   "TMe-2765" "TMe-1795" "TMe-2268" "TMe-3308"
#> [13] "TMe-1827" "TMe-1504" "TMe-3209" "TMe-2715" "TMe-2056" "TMe-369" 
#> [19] "TMe-74"   "TMe-3564" "TMe-1831" "TMe-3141" "TMe-2"    "TMe-2412"
#> [25] "TMe-509"  "TMe-3547" "TMe-47"   "TMe-3052" "TMe-3371" "TMe-3805"
#> [31] "TMe-160" 
#> 
#> [[3]]
#>  [1] "TMe-3143" "TMe-116"  "TMe-434"  "TMe-3299" "TMe-1863" "TMe-203" 
#>  [7] "TMe-3028" "TMe-1897" "TMe-3147" "TMe-635"  "TMe-3631" "TMe-2756"
#> [13] "TMe-3731" "TMe-174"  "TMe-1965" "TMe-2926" "TMe-1889" "TMe-3128"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-1162" "TMe-608"  "TMe-1020" "TMe-1577"
#>  [7] "TMe-209"  "TMe-3168" "TMe-875"  "TMe-3434" "TMe-3054" "TMe-368" 
#> [13] "TMe-87"   "TMe-3409" "TMe-1106" "TMe-3250" "TMe-3760" "TMe-3581"
#> [19] "TMe-82"   "TMe-1297" "TMe-380"  "TMe-2841" "TMe-3417" "TMe-2971"
#> [25] "TMe-215"  "TMe-3730" "TMe-3233" "TMe-210"  "TMe-2399" "TMe-43"  
#> [31] "TMe-2240" "TMe-2552" "TMe-3189" "TMe-1350" "TMe-1652"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1099" "TMe-167"  "TMe-376"  "TMe-2688"
#>  [7] "TMe-786"  "TMe-1873" "TMe-3411" "TMe-1551" "TMe-2192" "TMe-137" 
#> [13] "TMe-378"  "TMe-825"  "TMe-1366" "TMe-603"  "TMe-2905" "TMe-2753"
#> [19] "TMe-2057" "TMe-212"  "TMe-1382" "TMe-777"  "TMe-2400" "TMe-189" 
#> [25] "TMe-347"  "TMe-1455" "TMe-945"  "TMe-2213" "TMe-2907" "TMe-2853"
#> [31] "TMe-456"  "TMe-969"  "TMe-2003" "TMe-1290" "TMe-861"  "TMe-1924"
#> [37] "TMe-1220" "TMe-194"  "TMe-929"  "TMe-3659"
#> 
#> [[6]]
#>  [1] "TMe-1548" "TMe-1460" "TMe-1035" "TMe-705"  "TMe-1985" "TMe-1403"
#>  [7] "TMe-903"  "TMe-893"  "TMe-1554" "TMe-580"  "TMe-531"  "TMe-1518"
#> [13] "TMe-963"  "TMe-268"  "TMe-1232" "TMe-846"  "TMe-2035"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Pooled Brillouin diversity index")



# Mean Margalef's richness index
randomsel_mean_margalef <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_margalef,
                   metric = "mean", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_mean_margalef
#> [[1]]
#>  [1] "TMe-3423" "TMe-1964" "TMe-2255" "TMe-3398" "TMe-2912" "TMe-2083"
#>  [7] "TMe-3202" "TMe-3601" "TMe-3623" "TMe-1360" "TMe-839"  "TMe-2084"
#> [13] "TMe-3499" "TMe-3292" "TMe-3002" "TMe-1333" "TMe-500"  "TMe-2010"
#> [19] "TMe-3169" "TMe-3359" "TMe-815"  "TMe-888"  "TMe-2965" "TMe-3488"
#> [25] "TMe-2453" "TMe-1914" "TMe-1823" "TMe-670"  "TMe-3065" "TMe-566" 
#> [31] "TMe-3462"
#> 
#> [[2]]
#>  [1] "TMe-2757" "TMe-1323" "TMe-3805" "TMe-2995" "TMe-3466" "TMe-2998"
#>  [7] "TMe-3200" "TMe-1860" "TMe-1827" "TMe-2058" "TMe-2611" "TMe-3447"
#> [13] "TMe-2862" "TMe-589"  "TMe-1444" "TMe-2217" "TMe-409"  "TMe-3533"
#> [19] "TMe-681"  "TMe-74"   "TMe-2033" "TMe-2915" "TMe-2860" "TMe-2304"
#> [25] "TMe-3093" "TMe-289"  "TMe-3101" "TMe-3605" "TMe-601"  "TMe-1919"
#> [31] "TMe-1969"
#> 
#> [[3]]
#>  [1] "TMe-3772" "TMe-2481" "TMe-3631" "TMe-1868" "TMe-3592" "TMe-2285"
#>  [7] "TMe-3128" "TMe-2811" "TMe-2977" "TMe-1897" "TMe-2926" "TMe-3346"
#> [13] "TMe-3287" "TMe-3010" "TMe-1790" "TMe-3638" "TMe-64"   "TMe-3407"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-897"  "TMe-1123" "TMe-3256" "TMe-1155"
#>  [7] "TMe-3269" "TMe-3428" "TMe-108"  "TMe-5"    "TMe-3558" "TMe-467" 
#> [13] "TMe-318"  "TMe-596"  "TMe-1144" "TMe-3594" "TMe-1456" "TMe-317" 
#> [19] "TMe-3231" "TMe-2788" "TMe-3316" "TMe-3267" "TMe-556"  "TMe-279" 
#> [25] "TMe-427"  "TMe-2800" "TMe-2809" "TMe-394"  "TMe-885"  "TMe-1806"
#> [31] "TMe-3780" "TMe-3255" "TMe-2399" "TMe-2350" "TMe-1867"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-457"  "TMe-1311" "TMe-1440" "TMe-3628"
#>  [7] "TMe-1098" "TMe-1609" "TMe-700"  "TMe-2643" "TMe-1290" "TMe-2119"
#> [13] "TMe-712"  "TMe-3034" "TMe-326"  "TMe-863"  "TMe-2055" "TMe-2290"
#> [19] "TMe-2435" "TMe-651"  "TMe-1227" "TMe-2531" "TMe-2016" "TMe-344" 
#> [25] "TMe-347"  "TMe-1381" "TMe-1427" "TMe-1446" "TMe-585"  "TMe-835" 
#> [31] "TMe-1257" "TMe-412"  "TMe-355"  "TMe-669"  "TMe-1622" "TMe-2853"
#> [37] "TMe-481"  "TMe-397"  "TMe-1521" "TMe-378" 
#> 
#> [[6]]
#>  [1] "TMe-213"  "TMe-2510" "TMe-580"  "TMe-1232" "TMe-130"  "TMe-2983"
#>  [7] "TMe-1858" "TMe-985"  "TMe-845"  "TMe-690"  "TMe-3368" "TMe-683" 
#> [13] "TMe-1614" "TMe-1503" "TMe-188"  "TMe-791"  "TMe-3095"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Mean Margalef's diversity index")


# Pooled Margalef's richness index
randomsel_sum_margalef <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_margalef,
                   metric = "pooled", search = "random", local.search = NULL,
                   n.iter = 100)
randomsel_sum_margalef
#> [[1]]
#>  [1] "TMe-3423" "TMe-1191" "TMe-1247" "TMe-1869" "TMe-3236" "TMe-2010"
#>  [7] "TMe-1468" "TMe-3485" "TMe-952"  "TMe-2909" "TMe-3111" "TMe-3415"
#> [13] "TMe-3002" "TMe-1221" "TMe-3220" "TMe-3264" "TMe-489"  "TMe-20"  
#> [19] "TMe-2965" "TMe-3539" "TMe-2191" "TMe-2031" "TMe-3325" "TMe-2004"
#> [25] "TMe-694"  "TMe-2027" "TMe-1226" "TMe-3031" "TMe-2949" "TMe-1140"
#> [31] "TMe-3514"
#> 
#> [[2]]
#>  [1] "TMe-3466" "TMe-2903" "TMe-3140" "TMe-2757" "TMe-2268" "TMe-3498"
#>  [7] "TMe-3671" "TMe-1251" "TMe-2021" "TMe-1323" "TMe-2304" "TMe-1860"
#> [13] "TMe-3401" "TMe-2"    "TMe-369"  "TMe-2860" "TMe-2957" "TMe-1184"
#> [19] "TMe-404"  "TMe-2352" "TMe-2952" "TMe-2915" "TMe-2935" "TMe-1385"
#> [25] "TMe-2765" "TMe-1623" "TMe-821"  "TMe-3564" "TMe-3365" "TMe-2242"
#> [31] "TMe-3237"
#> 
#> [[3]]
#>  [1] "TMe-3575" "TMe-3731" "TMe-2926" "TMe-3592" "TMe-3046" "TMe-267" 
#>  [7] "TMe-1889" "TMe-2977" "TMe-70"   "TMe-1868" "TMe-1897" "TMe-3663"
#> [13] "TMe-3118" "TMe-2356" "TMe-1593" "TMe-3638" "TMe-2532" "TMe-104" 
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-3451" "TMe-2924" "TMe-279"  "TMe-241" 
#>  [7] "TMe-2776" "TMe-153"  "TMe-63"   "TMe-2367" "TMe-1128" "TMe-3248"
#> [13] "TMe-3511" "TMe-3527" "TMe-2036" "TMe-885"  "TMe-1480" "TMe-1350"
#> [19] "TMe-1328" "TMe-1923" "TMe-962"  "TMe-3581" "TMe-215"  "TMe-698" 
#> [25] "TMe-184"  "TMe-372"  "TMe-714"  "TMe-3327" "TMe-2979" "TMe-3073"
#> [31] "TMe-3312" "TMe-1948" "TMe-280"  "TMe-1525" "TMe-1700"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-747"  "TMe-950"  "TMe-532"  "TMe-1307"
#>  [7] "TMe-475"  "TMe-2589" "TMe-412"  "TMe-723"  "TMe-969"  "TMe-1098"
#> [13] "TMe-193"  "TMe-1762" "TMe-795"  "TMe-2542" "TMe-1381" "TMe-665" 
#> [19] "TMe-510"  "TMe-2015" "TMe-197"  "TMe-439"  "TMe-1367" "TMe-2820"
#> [25] "TMe-1311" "TMe-168"  "TMe-464"  "TMe-2041" "TMe-1011" "TMe-2326"
#> [31] "TMe-948"  "TMe-729"  "TMe-336"  "TMe-487"  "TMe-256"  "TMe-1812"
#> [37] "TMe-2119" "TMe-3185" "TMe-1391" "TMe-1482"
#> 
#> [[6]]
#>  [1] "TMe-625"  "TMe-1723" "TMe-856"  "TMe-1035" "TMe-903"  "TMe-1124"
#>  [7] "TMe-1661" "TMe-598"  "TMe-626"  "TMe-2983" "TMe-315"  "TMe-725" 
#> [13] "TMe-3095" "TMe-838"  "TMe-693"  "TMe-652"  "TMe-1775"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Random search",
       subtitle = "Pooled Margalef's diversity index")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Greedy search with 1-opt best improvement
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mean richness
greedysel_best_mean_richness <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "richness",
                   metric = "mean", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_mean_richness
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-677"  "TMe-3475" "TMe-1717"
#>  [7] "TMe-841"  "TMe-1425" "TMe-2010" "TMe-3030" "TMe-3292" "TMe-2069"
#> [13] "TMe-132"  "TMe-778"  "TMe-1226" "TMe-1360" "TMe-1922" "TMe-3262"
#> [19] "TMe-1915" "TMe-41"   "TMe-44"   "TMe-117"  "TMe-290"  "TMe-299" 
#> [25] "TMe-300"  "TMe-407"  "TMe-410"  "TMe-438"  "TMe-446"  "TMe-469" 
#> [31] "TMe-489" 
#> 
#> [[2]]
#>  [1] "TMe-1732" "TMe-369"  "TMe-1919" "TMe-1251" "TMe-3366" "TMe-1907"
#>  [7] "TMe-160"  "TMe-2033" "TMe-3200" "TMe-289"  "TMe-509"  "TMe-1241"
#> [13] "TMe-2304" "TMe-3466" "TMe-47"   "TMe-1461" "TMe-3222" "TMe-3237"
#> [19] "TMe-2"    "TMe-39"   "TMe-60"   "TMe-85"   "TMe-339"  "TMe-365" 
#> [25] "TMe-377"  "TMe-404"  "TMe-409"  "TMe-433"  "TMe-450"  "TMe-455" 
#> [31] "TMe-477" 
#> 
#> [[3]]
#>  [1] "TMe-64"   "TMe-1868" "TMe-785"  "TMe-946"  "TMe-1993" "TMe-3638"
#>  [7] "TMe-1897" "TMe-1819" "TMe-267"  "TMe-4"    "TMe-13"   "TMe-123" 
#> [13] "TMe-2592" "TMe-3046" "TMe-35"   "TMe-2926" "TMe-3596" "TMe-6"   
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1350" "TMe-552"  "TMe-1434"
#>  [7] "TMe-3255" "TMe-3527" "TMe-3054" "TMe-540"  "TMe-3428" "TMe-241" 
#> [13] "TMe-103"  "TMe-266"  "TMe-427"  "TMe-460"  "TMe-550"  "TMe-2809"
#> [19] "TMe-3273" "TMe-1020" "TMe-1511" "TMe-2567" "TMe-11"   "TMe-12"  
#> [25] "TMe-57"   "TMe-62"   "TMe-90"   "TMe-107"  "TMe-108"  "TMe-150" 
#> [31] "TMe-154"  "TMe-170"  "TMe-182"  "TMe-191"  "TMe-204" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2307" "TMe-3185" "TMe-197"  "TMe-723" 
#>  [7] "TMe-378"  "TMe-2953" "TMe-277"  "TMe-2853" "TMe-1004" "TMe-336" 
#> [13] "TMe-399"  "TMe-457"  "TMe-1904" "TMe-2041" "TMe-181"  "TMe-48"  
#> [19] "TMe-98"   "TMe-99"   "TMe-158"  "TMe-167"  "TMe-212"  "TMe-275" 
#> [25] "TMe-294"  "TMe-305"  "TMe-307"  "TMe-332"  "TMe-334"  "TMe-340" 
#> [31] "TMe-347"  "TMe-359"  "TMe-360"  "TMe-362"  "TMe-363"  "TMe-373" 
#> [37] "TMe-376"  "TMe-382"  "TMe-385"  "TMe-390" 
#> 
#> [[6]]
#>  [1] "TMe-1633" "TMe-2983" "TMe-1171" "TMe-2510" "TMe-315"  "TMe-2534"
#>  [7] "TMe-1445" "TMe-130"  "TMe-1875" "TMe-693"  "TMe-751"  "TMe-1174"
#> [13] "TMe-1353" "TMe-1775" "TMe-856"  "TMe-3095" "TMe-138" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Mean richness")


# Pooled richness
greedysel_best_sum_richness <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "richness",
                   metric = "pooled", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_sum_richness
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-677"  "TMe-3475" "TMe-1717"
#>  [7] "TMe-841"  "TMe-1425" "TMe-2010" "TMe-3030" "TMe-3292" "TMe-2069"
#> [13] "TMe-132"  "TMe-778"  "TMe-1226" "TMe-1360" "TMe-1922" "TMe-3262"
#> [19] "TMe-1915" "TMe-41"   "TMe-44"   "TMe-117"  "TMe-290"  "TMe-299" 
#> [25] "TMe-300"  "TMe-407"  "TMe-410"  "TMe-438"  "TMe-446"  "TMe-469" 
#> [31] "TMe-489" 
#> 
#> [[2]]
#>  [1] "TMe-3222" "TMe-47"   "TMe-1827" "TMe-2903" "TMe-509"  "TMe-477" 
#>  [7] "TMe-3200" "TMe-289"  "TMe-2033" "TMe-1907" "TMe-160"  "TMe-369" 
#> [13] "TMe-1241" "TMe-1385" "TMe-2304" "TMe-3466" "TMe-1461" "TMe-3237"
#> [19] "TMe-2"    "TMe-39"   "TMe-60"   "TMe-85"   "TMe-339"  "TMe-365" 
#> [25] "TMe-377"  "TMe-404"  "TMe-409"  "TMe-433"  "TMe-450"  "TMe-455" 
#> [31] "TMe-528" 
#> 
#> [[3]]
#>  [1] "TMe-3335" "TMe-3596" "TMe-785"  "TMe-946"  "TMe-261"  "TMe-13"  
#>  [7] "TMe-3048" "TMe-420"  "TMe-1868" "TMe-1993" "TMe-4"    "TMe-123" 
#> [13] "TMe-267"  "TMe-2592" "TMe-3046" "TMe-35"   "TMe-2926" "TMe-6"   
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1350" "TMe-552"  "TMe-1434"
#>  [7] "TMe-3255" "TMe-3527" "TMe-3054" "TMe-540"  "TMe-3428" "TMe-241" 
#> [13] "TMe-103"  "TMe-266"  "TMe-427"  "TMe-460"  "TMe-550"  "TMe-2809"
#> [19] "TMe-3273" "TMe-1020" "TMe-1511" "TMe-2567" "TMe-11"   "TMe-12"  
#> [25] "TMe-57"   "TMe-62"   "TMe-90"   "TMe-107"  "TMe-108"  "TMe-150" 
#> [31] "TMe-154"  "TMe-170"  "TMe-182"  "TMe-191"  "TMe-204" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2307" "TMe-3185" "TMe-197"  "TMe-723" 
#>  [7] "TMe-378"  "TMe-2953" "TMe-277"  "TMe-2853" "TMe-1004" "TMe-336" 
#> [13] "TMe-399"  "TMe-457"  "TMe-1904" "TMe-2041" "TMe-181"  "TMe-48"  
#> [19] "TMe-98"   "TMe-99"   "TMe-158"  "TMe-167"  "TMe-212"  "TMe-275" 
#> [25] "TMe-294"  "TMe-305"  "TMe-307"  "TMe-332"  "TMe-334"  "TMe-340" 
#> [31] "TMe-347"  "TMe-359"  "TMe-360"  "TMe-362"  "TMe-363"  "TMe-373" 
#> [37] "TMe-376"  "TMe-382"  "TMe-385"  "TMe-390" 
#> 
#> [[6]]
#>  [1] "TMe-693"  "TMe-514"  "TMe-3549" "TMe-1079" "TMe-315"  "TMe-2510"
#>  [7] "TMe-1883" "TMe-2064" "TMe-213"  "TMe-1174" "TMe-1392" "TMe-130" 
#> [13] "TMe-1775" "TMe-1875" "TMe-856"  "TMe-3095" "TMe-138" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Pooled richness")


# Mean Shannon-Weaver diversity index
greedysel_best_mean_shannon <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "shannon",
                   metric = "mean", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_mean_shannon
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-3394" "TMe-2069" "TMe-1425"
#>  [7] "TMe-41"   "TMe-2010" "TMe-3292" "TMe-566"  "TMe-3030" "TMe-1360"
#> [13] "TMe-2993" "TMe-2965" "TMe-1960" "TMe-2967" "TMe-3319" "TMe-1423"
#> [19] "TMe-606"  "TMe-1190" "TMe-3534" "TMe-3396" "TMe-2785" "TMe-2237"
#> [25] "TMe-3112" "TMe-1922" "TMe-3102" "TMe-1096" "TMe-778"  "TMe-1786"
#> [31] "TMe-3475"
#> 
#> [[2]]
#>  [1] "TMe-3222" "TMe-369"  "TMe-2033" "TMe-2952" "TMe-2056" "TMe-47"  
#>  [7] "TMe-509"  "TMe-3690" "TMe-2866" "TMe-3338" "TMe-3141" "TMe-3200"
#> [13] "TMe-1919" "TMe-3771" "TMe-3466" "TMe-3237" "TMe-1907" "TMe-1754"
#> [19] "TMe-1444" "TMe-3286" "TMe-289"  "TMe-1107" "TMe-1461" "TMe-2128"
#> [25] "TMe-2860" "TMe-1385" "TMe-2329" "TMe-160"  "TMe-40"   "TMe-3805"
#> [31] "TMe-1241"
#> 
#> [[3]]
#>  [1] "TMe-234"  "TMe-1897" "TMe-2481" "TMe-2823" "TMe-785"  "TMe-1993"
#>  [7] "TMe-3596" "TMe-3274" "TMe-946"  "TMe-1868" "TMe-1593" "TMe-878" 
#> [13] "TMe-116"  "TMe-3133" "TMe-2926" "TMe-3383" "TMe-3445" "TMe-3715"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-812"  "TMe-1434" "TMe-3255"
#>  [7] "TMe-3422" "TMe-1020" "TMe-3730" "TMe-3167" "TMe-225"  "TMe-2809"
#> [13] "TMe-1988" "TMe-1511" "TMe-3527" "TMe-3442" "TMe-3054" "TMe-3225"
#> [19] "TMe-2971" "TMe-241"  "TMe-2924" "TMe-1376" "TMe-2956" "TMe-3273"
#> [25] "TMe-2788" "TMe-266"  "TMe-1525" "TMe-3242" "TMe-3542" "TMe-550" 
#> [31] "TMe-2567" "TMe-3760" "TMe-1123" "TMe-1419" "TMe-663" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-208"  "TMe-3185" "TMe-3411" "TMe-3034"
#>  [7] "TMe-277"  "TMe-2953" "TMe-378"  "TMe-723"  "TMe-2041" "TMe-2853"
#> [13] "TMe-2589" "TMe-2531" "TMe-284"  "TMe-1381" "TMe-1458" "TMe-969" 
#> [19] "TMe-712"  "TMe-1004" "TMe-2121" "TMe-2355" "TMe-2790" "TMe-287" 
#> [25] "TMe-336"  "TMe-3356" "TMe-1500" "TMe-536"  "TMe-2820" "TMe-3332"
#> [31] "TMe-1963" "TMe-294"  "TMe-3628" "TMe-759"  "TMe-2612" "TMe-3311"
#> [37] "TMe-1257" "TMe-2863" "TMe-1099" "TMe-2907"
#> 
#> [[6]]
#>  [1] "TMe-751"  "TMe-2064" "TMe-1509" "TMe-1646" "TMe-3549" "TMe-693" 
#>  [7] "TMe-2067" "TMe-1875" "TMe-1232" "TMe-1428" "TMe-1174" "TMe-1566"
#> [13] "TMe-130"  "TMe-2035" "TMe-3095" "TMe-1858" "TMe-580" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Mean Shannon-Weaver diversity index")


# Pooled Shannon-Weaver diversity index
greedysel_best_sum_shannon <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "shannon",
                   metric = "pooled", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_sum_shannon
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-3394" "TMe-2069" "TMe-1425"
#>  [7] "TMe-41"   "TMe-2010" "TMe-3292" "TMe-566"  "TMe-3030" "TMe-1360"
#> [13] "TMe-2993" "TMe-2965" "TMe-1960" "TMe-2967" "TMe-3319" "TMe-1423"
#> [19] "TMe-606"  "TMe-1190" "TMe-3534" "TMe-3396" "TMe-2785" "TMe-2237"
#> [25] "TMe-3112" "TMe-1922" "TMe-3102" "TMe-1096" "TMe-778"  "TMe-1786"
#> [31] "TMe-3475"
#> 
#> [[2]]
#>  [1] "TMe-3766" "TMe-369"  "TMe-2128" "TMe-3200" "TMe-2033" "TMe-2952"
#>  [7] "TMe-3771" "TMe-1907" "TMe-3466" "TMe-3338" "TMe-2866" "TMe-47"  
#> [13] "TMe-509"  "TMe-2056" "TMe-2329" "TMe-3141" "TMe-40"   "TMe-3286"
#> [19] "TMe-289"  "TMe-1919" "TMe-1444" "TMe-1754" "TMe-1385" "TMe-2860"
#> [25] "TMe-1461" "TMe-1107" "TMe-3690" "TMe-160"  "TMe-1241" "TMe-3009"
#> [31] "TMe-3237"
#> 
#> [[3]]
#>  [1] "TMe-35"   "TMe-200"  "TMe-1993" "TMe-1897" "TMe-3383" "TMe-116" 
#>  [7] "TMe-878"  "TMe-2926" "TMe-3565" "TMe-2481" "TMe-3596" "TMe-3274"
#> [13] "TMe-3445" "TMe-234"  "TMe-785"  "TMe-2532" "TMe-2823" "TMe-1868"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-812"  "TMe-1434" "TMe-3255"
#>  [7] "TMe-3422" "TMe-1020" "TMe-3730" "TMe-3167" "TMe-225"  "TMe-2809"
#> [13] "TMe-1988" "TMe-1511" "TMe-3527" "TMe-3442" "TMe-3054" "TMe-3225"
#> [19] "TMe-2971" "TMe-241"  "TMe-2924" "TMe-1376" "TMe-2956" "TMe-3273"
#> [25] "TMe-2788" "TMe-266"  "TMe-1525" "TMe-3242" "TMe-3542" "TMe-550" 
#> [31] "TMe-2567" "TMe-3760" "TMe-1123" "TMe-1419" "TMe-663" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-208"  "TMe-3185" "TMe-3411" "TMe-3034"
#>  [7] "TMe-277"  "TMe-2953" "TMe-378"  "TMe-723"  "TMe-2041" "TMe-2853"
#> [13] "TMe-2589" "TMe-2531" "TMe-284"  "TMe-1381" "TMe-1458" "TMe-969" 
#> [19] "TMe-712"  "TMe-1004" "TMe-2121" "TMe-2355" "TMe-2790" "TMe-287" 
#> [25] "TMe-336"  "TMe-3356" "TMe-1500" "TMe-536"  "TMe-2820" "TMe-3332"
#> [31] "TMe-1963" "TMe-294"  "TMe-3628" "TMe-759"  "TMe-2612" "TMe-3311"
#> [37] "TMe-1257" "TMe-2863" "TMe-1099" "TMe-2907"
#> 
#> [[6]]
#>  [1] "TMe-1883" "TMe-3095" "TMe-1124" "TMe-2791" "TMe-1428" "TMe-1232"
#>  [7] "TMe-1566" "TMe-2035" "TMe-1509" "TMe-1174" "TMe-3549" "TMe-1875"
#> [13] "TMe-1858" "TMe-580"  "TMe-751"  "TMe-2067" "TMe-693" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Pooled Shannon-Weaver diversity index")


# Mean Gini-Simpson diversity index
greedysel_best_mean_simpson <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "simpson",
                   metric = "mean", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_mean_simpson
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-3452" "TMe-2069" "TMe-1425"
#>  [7] "TMe-3389" "TMe-1890" "TMe-41"   "TMe-3292" "TMe-1360" "TMe-1960"
#> [13] "TMe-2965" "TMe-1190" "TMe-2967" "TMe-3030" "TMe-1423" "TMe-778" 
#> [19] "TMe-606"  "TMe-3475" "TMe-2993" "TMe-3102" "TMe-2372" "TMe-841" 
#> [25] "TMe-3319" "TMe-132"  "TMe-3208" "TMe-1086" "TMe-3705" "TMe-1922"
#> [31] "TMe-1564"
#> 
#> [[2]]
#>  [1] "TMe-3286" "TMe-3141" "TMe-1919" "TMe-3009" "TMe-3338" "TMe-3466"
#>  [7] "TMe-369"  "TMe-2033" "TMe-3200" "TMe-2866" "TMe-47"   "TMe-1754"
#> [13] "TMe-1907" "TMe-1385" "TMe-3771" "TMe-2952" "TMe-3766" "TMe-2056"
#> [19] "TMe-509"  "TMe-1444" "TMe-1107" "TMe-86"   "TMe-3237" "TMe-2128"
#> [25] "TMe-3805" "TMe-3690" "TMe-478"  "TMe-160"  "TMe-40"   "TMe-2329"
#> [31] "TMe-2860"
#> 
#> [[3]]
#>  [1] "TMe-3804" "TMe-381"  "TMe-1868" "TMe-3596" "TMe-1897" "TMe-3383"
#>  [7] "TMe-2977" "TMe-267"  "TMe-785"  "TMe-2926" "TMe-200"  "TMe-1593"
#> [13] "TMe-3445" "TMe-3274" "TMe-3128" "TMe-1200" "TMe-141"  "TMe-1993"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1434" "TMe-388"  "TMe-3275"
#>  [7] "TMe-1083" "TMe-3225" "TMe-2809" "TMe-2924" "TMe-1020" "TMe-962" 
#> [13] "TMe-812"  "TMe-1511" "TMe-3442" "TMe-3054" "TMe-2971" "TMe-1988"
#> [19] "TMe-241"  "TMe-2552" "TMe-3417" "TMe-3542" "TMe-3527" "TMe-1525"
#> [25] "TMe-650"  "TMe-180"  "TMe-2240" "TMe-1419" "TMe-550"  "TMe-3242"
#> [31] "TMe-3760" "TMe-3198" "TMe-2956" "TMe-3273" "TMe-3390"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1257" "TMe-3185" "TMe-712"  "TMe-1381"
#>  [7] "TMe-2589" "TMe-1963" "TMe-723"  "TMe-1458" "TMe-2953" "TMe-2531"
#> [13] "TMe-378"  "TMe-284"  "TMe-1534" "TMe-3034" "TMe-294"  "TMe-969" 
#> [19] "TMe-2355" "TMe-1500" "TMe-3356" "TMe-2820" "TMe-536"  "TMe-2041"
#> [25] "TMe-759"  "TMe-2003" "TMe-2853" "TMe-2121" "TMe-3311" "TMe-1099"
#> [31] "TMe-1357" "TMe-2612" "TMe-3332" "TMe-336"  "TMe-158"  "TMe-2907"
#> [37] "TMe-2790" "TMe-3628" "TMe-48"   "TMe-277" 
#> 
#> [[6]]
#>  [1] "TMe-1428" "TMe-2510" "TMe-2064" "TMe-1353" "TMe-751"  "TMe-1232"
#>  [7] "TMe-3095" "TMe-1035" "TMe-361"  "TMe-236"  "TMe-693"  "TMe-580" 
#> [13] "TMe-3549" "TMe-1875" "TMe-2791" "TMe-470"  "TMe-2983"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Mean Gini-Simpson diversity index")


# Pooled Gini-Simpson diversity index
greedysel_best_sum_simpson <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "simpson",
                   metric = "pooled", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_sum_simpson
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-3282" "TMe-41"   "TMe-839"  "TMe-1960"
#>  [7] "TMe-2965" "TMe-2967" "TMe-1423" "TMe-1190" "TMe-3292" "TMe-3437"
#> [13] "TMe-1425" "TMe-2785" "TMe-867"  "TMe-3394" "TMe-3102" "TMe-2993"
#> [19] "TMe-1360" "TMe-3475" "TMe-841"  "TMe-2462" "TMe-28"   "TMe-2069"
#> [25] "TMe-606"  "TMe-717"  "TMe-1600" "TMe-3151" "TMe-3319" "TMe-1786"
#> [31] "TMe-1333"
#> 
#> [[2]]
#>  [1] "TMe-3690" "TMe-1444" "TMe-2866" "TMe-47"   "TMe-3141" "TMe-3338"
#>  [7] "TMe-2056" "TMe-3009" "TMe-3466" "TMe-369"  "TMe-3237" "TMe-1907"
#> [13] "TMe-1107" "TMe-509"  "TMe-3200" "TMe-2033" "TMe-1754" "TMe-2952"
#> [19] "TMe-2128" "TMe-3771" "TMe-86"   "TMe-1385" "TMe-1919" "TMe-3766"
#> [25] "TMe-40"   "TMe-3805" "TMe-478"  "TMe-160"  "TMe-3286" "TMe-2329"
#> [31] "TMe-2860"
#> 
#> [[3]]
#>  [1] "TMe-878"  "TMe-2926" "TMe-3383" "TMe-200"  "TMe-3596" "TMe-1897"
#>  [7] "TMe-3274" "TMe-1593" "TMe-3565" "TMe-2481" "TMe-267"  "TMe-3445"
#> [13] "TMe-773"  "TMe-3715" "TMe-116"  "TMe-141"  "TMe-1868" "TMe-785" 
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1434" "TMe-388"  "TMe-3275"
#>  [7] "TMe-1083" "TMe-3225" "TMe-2809" "TMe-2924" "TMe-1020" "TMe-1511"
#> [13] "TMe-3707" "TMe-2971" "TMe-3442" "TMe-550"  "TMe-962"  "TMe-241" 
#> [19] "TMe-2956" "TMe-2552" "TMe-1406" "TMe-3730" "TMe-3527" "TMe-1525"
#> [25] "TMe-3542" "TMe-3054" "TMe-650"  "TMe-1988" "TMe-3167" "TMe-3417"
#> [31] "TMe-180"  "TMe-3273" "TMe-3390" "TMe-3576" "TMe-225" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2435" "TMe-3185" "TMe-712"  "TMe-1381"
#>  [7] "TMe-2589" "TMe-2750" "TMe-723"  "TMe-1458" "TMe-2953" "TMe-2531"
#> [13] "TMe-378"  "TMe-284"  "TMe-2863" "TMe-3034" "TMe-1357" "TMe-2820"
#> [19] "TMe-1500" "TMe-3356" "TMe-294"  "TMe-3311" "TMe-2853" "TMe-2121"
#> [25] "TMe-2041" "TMe-158"  "TMe-536"  "TMe-969"  "TMe-759"  "TMe-208" 
#> [31] "TMe-2003" "TMe-1099" "TMe-2907" "TMe-336"  "TMe-3628" "TMe-2612"
#> [37] "TMe-48"   "TMe-2790" "TMe-2355" "TMe-277" 
#> 
#> [[6]]
#>  [1] "TMe-1403" "TMe-693"  "TMe-1509" "TMe-3549" "TMe-1079" "TMe-3095"
#>  [7] "TMe-1353" "TMe-2064" "TMe-1232" "TMe-1053" "TMe-1035" "TMe-2791"
#> [13] "TMe-1875" "TMe-1678" "TMe-2510" "TMe-1124" "TMe-2067"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Pooled Gini-Simpson diversity index")


# Mean McIntosh diversity index
greedysel_best_mean_mcintosh <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "mcintosh",
                   metric = "pooled", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_mean_mcintosh
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-3208" "TMe-1468" "TMe-2967"
#>  [7] "TMe-1425" "TMe-2785" "TMe-1190" "TMe-3437" "TMe-1086" "TMe-2934"
#> [13] "TMe-2993" "TMe-1360" "TMe-2965" "TMe-3433" "TMe-407"  "TMe-3292"
#> [19] "TMe-41"   "TMe-1960" "TMe-2069" "TMe-1423" "TMe-3102" "TMe-606" 
#> [25] "TMe-3396" "TMe-28"   "TMe-3236" "TMe-1333" "TMe-3319" "TMe-3030"
#> [31] "TMe-3705"
#> 
#> [[2]]
#>  [1] "TMe-478"  "TMe-3140" "TMe-47"   "TMe-1385" "TMe-2128" "TMe-3141"
#>  [7] "TMe-3338" "TMe-3466" "TMe-369"  "TMe-1754" "TMe-1107" "TMe-1907"
#> [13] "TMe-160"  "TMe-40"   "TMe-509"  "TMe-1919" "TMe-3805" "TMe-3200"
#> [19] "TMe-1444" "TMe-2056" "TMe-2952" "TMe-2033" "TMe-2866" "TMe-3771"
#> [25] "TMe-3766" "TMe-2329" "TMe-2860" "TMe-3286" "TMe-636"  "TMe-3690"
#> [31] "TMe-3237"
#> 
#> [[3]]
#>  [1] "TMe-2897" "TMe-200"  "TMe-2926" "TMe-785"  "TMe-3383" "TMe-1897"
#>  [7] "TMe-3715" "TMe-2481" "TMe-3596" "TMe-3565" "TMe-3274" "TMe-773" 
#> [13] "TMe-2823" "TMe-3445" "TMe-116"  "TMe-1868" "TMe-1863" "TMe-141" 
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-2552" "TMe-1434" "TMe-3255"
#>  [7] "TMe-3442" "TMe-2971" "TMe-209"  "TMe-2809" "TMe-1575" "TMe-1083"
#> [13] "TMe-2924" "TMe-3707" "TMe-1511" "TMe-241"  "TMe-388"  "TMe-2956"
#> [19] "TMe-550"  "TMe-3225" "TMe-650"  "TMe-3422" "TMe-180"  "TMe-225" 
#> [25] "TMe-1020" "TMe-1597" "TMe-3054" "TMe-3527" "TMe-3542" "TMe-1525"
#> [31] "TMe-1376" "TMe-3417" "TMe-3390" "TMe-1988" "TMe-2240"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2790" "TMe-3185" "TMe-277"  "TMe-1381"
#>  [7] "TMe-378"  "TMe-2953" "TMe-723"  "TMe-2589" "TMe-2531" "TMe-344" 
#> [13] "TMe-2121" "TMe-2853" "TMe-712"  "TMe-3332" "TMe-1924" "TMe-3034"
#> [19] "TMe-1500" "TMe-1357" "TMe-208"  "TMe-969"  "TMe-287"  "TMe-284" 
#> [25] "TMe-3356" "TMe-294"  "TMe-759"  "TMe-536"  "TMe-2355" "TMe-2820"
#> [31] "TMe-2612" "TMe-2041" "TMe-1963" "TMe-3628" "TMe-1458" "TMe-3311"
#> [37] "TMe-2003" "TMe-1099" "TMe-2863" "TMe-755" 
#> 
#> [[6]]
#>  [1] "TMe-725"  "TMe-1403" "TMe-580"  "TMe-1124" "TMe-2510" "TMe-1353"
#>  [7] "TMe-3095" "TMe-1678" "TMe-2064" "TMe-3549" "TMe-693"  "TMe-1509"
#> [13] "TMe-1875" "TMe-2791" "TMe-2067" "TMe-1035" "TMe-1232"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Mean McIntosh diversity index")


# Pooled McIntosh diversity index
greedysel_best_sum_mcintosh <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "mcintosh",
                   metric = "pooled", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_sum_mcintosh
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-3208" "TMe-1468" "TMe-2967"
#>  [7] "TMe-1425" "TMe-2785" "TMe-1190" "TMe-3437" "TMe-1086" "TMe-2934"
#> [13] "TMe-2993" "TMe-1360" "TMe-2965" "TMe-3433" "TMe-407"  "TMe-3292"
#> [19] "TMe-41"   "TMe-1960" "TMe-2069" "TMe-1423" "TMe-3102" "TMe-606" 
#> [25] "TMe-3396" "TMe-28"   "TMe-3236" "TMe-1333" "TMe-3319" "TMe-3030"
#> [31] "TMe-3705"
#> 
#> [[2]]
#>  [1] "TMe-47"   "TMe-3009" "TMe-3805" "TMe-3466" "TMe-509"  "TMe-1919"
#>  [7] "TMe-369"  "TMe-1754" "TMe-3690" "TMe-3338" "TMe-3141" "TMe-2866"
#> [13] "TMe-2033" "TMe-3200" "TMe-2056" "TMe-1444" "TMe-1107" "TMe-86"  
#> [19] "TMe-3771" "TMe-2952" "TMe-3237" "TMe-1907" "TMe-2128" "TMe-1385"
#> [25] "TMe-160"  "TMe-40"   "TMe-2329" "TMe-2860" "TMe-478"  "TMe-3286"
#> [31] "TMe-3766"
#> 
#> [[3]]
#>  [1] "TMe-2532" "TMe-116"  "TMe-1897" "TMe-3638" "TMe-2823" "TMe-785" 
#>  [7] "TMe-3596" "TMe-946"  "TMe-3274" "TMe-3631" "TMe-234"  "TMe-1838"
#> [13] "TMe-3143" "TMe-2926" "TMe-3383" "TMe-1200" "TMe-200"  "TMe-2897"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-2552" "TMe-1434" "TMe-3255"
#>  [7] "TMe-3442" "TMe-2971" "TMe-209"  "TMe-2809" "TMe-1575" "TMe-1083"
#> [13] "TMe-2924" "TMe-3707" "TMe-1511" "TMe-241"  "TMe-388"  "TMe-2956"
#> [19] "TMe-550"  "TMe-3225" "TMe-650"  "TMe-3422" "TMe-180"  "TMe-225" 
#> [25] "TMe-1020" "TMe-1597" "TMe-3054" "TMe-3527" "TMe-3542" "TMe-1525"
#> [31] "TMe-1376" "TMe-3417" "TMe-3390" "TMe-1988" "TMe-2240"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2790" "TMe-3185" "TMe-277"  "TMe-1381"
#>  [7] "TMe-378"  "TMe-2953" "TMe-723"  "TMe-2589" "TMe-2531" "TMe-344" 
#> [13] "TMe-2121" "TMe-2853" "TMe-712"  "TMe-3332" "TMe-1924" "TMe-3034"
#> [19] "TMe-1500" "TMe-1357" "TMe-208"  "TMe-969"  "TMe-287"  "TMe-284" 
#> [25] "TMe-3356" "TMe-294"  "TMe-759"  "TMe-536"  "TMe-2355" "TMe-2820"
#> [31] "TMe-2612" "TMe-2041" "TMe-1963" "TMe-3628" "TMe-1458" "TMe-3311"
#> [37] "TMe-2003" "TMe-1099" "TMe-2863" "TMe-755" 
#> 
#> [[6]]
#>  [1] "TMe-1232" "TMe-2035" "TMe-1403" "TMe-470"  "TMe-2067" "TMe-693" 
#>  [7] "TMe-1428" "TMe-580"  "TMe-3095" "TMe-1035" "TMe-1875" "TMe-751" 
#> [13] "TMe-2510" "TMe-3549" "TMe-1353" "TMe-2791" "TMe-2064"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Pooled McIntosh diversity index")


# Mean Brillouin diversity index
greedysel_best_mean_brillouin <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_brillouin,
                   metric = "mean", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_mean_brillouin
#> [[1]]
#>  [1] "TMe-3423" "TMe-2084" "TMe-867"  "TMe-3475" "TMe-2069" "TMe-1425"
#>  [7] "TMe-3030" "TMe-1360" "TMe-3389" "TMe-3292" "TMe-41"   "TMe-2010"
#> [13] "TMe-1190" "TMe-2965" "TMe-1333" "TMe-2785" "TMe-2993" "TMe-1423"
#> [19] "TMe-2967" "TMe-606"  "TMe-2237" "TMe-1960" "TMe-3396" "TMe-1922"
#> [25] "TMe-3102" "TMe-778"  "TMe-132"  "TMe-2934" "TMe-1786" "TMe-3394"
#> [31] "TMe-566" 
#> 
#> [[2]]
#>  [1] "TMe-3237" "TMe-2866" "TMe-3766" "TMe-369"  "TMe-2033" "TMe-3338"
#>  [7] "TMe-3466" "TMe-1919" "TMe-1754" "TMe-1385" "TMe-3141" "TMe-47"  
#> [13] "TMe-1907" "TMe-3605" "TMe-509"  "TMe-40"   "TMe-2860" "TMe-1241"
#> [19] "TMe-3200" "TMe-2128" "TMe-3771" "TMe-3286" "TMe-2056" "TMe-2329"
#> [25] "TMe-160"  "TMe-2952" "TMe-289"  "TMe-1461" "TMe-1107" "TMe-3690"
#> [31] "TMe-1444"
#> 
#> [[3]]
#>  [1] "TMe-1868" "TMe-35"   "TMe-785"  "TMe-3383" "TMe-1897" "TMe-104" 
#>  [7] "TMe-2823" "TMe-2532" "TMe-1993" "TMe-3596" "TMe-3274" "TMe-3445"
#> [13] "TMe-878"  "TMe-2926" "TMe-1838" "TMe-234"  "TMe-200"  "TMe-1593"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1434" "TMe-3527" "TMe-3255"
#>  [7] "TMe-3542" "TMe-650"  "TMe-2971" "TMe-2809" "TMe-663"  "TMe-1419"
#> [13] "TMe-1511" "TMe-3054" "TMe-3442" "TMe-1988" "TMe-241"  "TMe-2924"
#> [19] "TMe-1123" "TMe-2567" "TMe-3225" "TMe-3730" "TMe-3242" "TMe-550" 
#> [25] "TMe-1020" "TMe-225"  "TMe-2240" "TMe-3422" "TMe-1525" "TMe-3273"
#> [31] "TMe-1376" "TMe-2956" "TMe-3167" "TMe-388"  "TMe-3760"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2907" "TMe-3185" "TMe-277"  "TMe-1381"
#>  [7] "TMe-378"  "TMe-2953" "TMe-723"  "TMe-2589" "TMe-2531" "TMe-1257"
#> [13] "TMe-2853" "TMe-2121" "TMe-2041" "TMe-969"  "TMe-3034" "TMe-712" 
#> [19] "TMe-287"  "TMe-2790" "TMe-1500" "TMe-3356" "TMe-2355" "TMe-284" 
#> [25] "TMe-1004" "TMe-48"   "TMe-2820" "TMe-336"  "TMe-294"  "TMe-3628"
#> [31] "TMe-158"  "TMe-536"  "TMe-1458" "TMe-759"  "TMe-2612" "TMe-1963"
#> [37] "TMe-3311" "TMe-2296" "TMe-208"  "TMe-1357"
#> 
#> [[6]]
#>  [1] "TMe-1566" "TMe-1079" "TMe-1509" "TMe-3095" "TMe-751"  "TMe-1174"
#>  [7] "TMe-1232" "TMe-1428" "TMe-2035" "TMe-130"  "TMe-693"  "TMe-1883"
#> [13] "TMe-2510" "TMe-3549" "TMe-1875" "TMe-315"  "TMe-2064"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Mean Brillouin diversity index")


# Pooled Brillouin diversity index
greedysel_best_sum_brillouin <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_brillouin,
                   metric = "pooled", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_sum_brillouin
#> [[1]]
#>  [1] "TMe-3423" "TMe-2084" "TMe-867"  "TMe-3475" "TMe-2069" "TMe-1425"
#>  [7] "TMe-3030" "TMe-1360" "TMe-3389" "TMe-3292" "TMe-41"   "TMe-2010"
#> [13] "TMe-1190" "TMe-2965" "TMe-1333" "TMe-2785" "TMe-2993" "TMe-1423"
#> [19] "TMe-2967" "TMe-606"  "TMe-2237" "TMe-1960" "TMe-3396" "TMe-1922"
#> [25] "TMe-3102" "TMe-778"  "TMe-132"  "TMe-2934" "TMe-1786" "TMe-3394"
#> [31] "TMe-566" 
#> 
#> [[2]]
#>  [1] "TMe-289"  "TMe-3805" "TMe-47"   "TMe-1461" "TMe-648"  "TMe-1919"
#>  [7] "TMe-3690" "TMe-1754" "TMe-3466" "TMe-369"  "TMe-2033" "TMe-3338"
#> [13] "TMe-3141" "TMe-509"  "TMe-2128" "TMe-160"  "TMe-40"   "TMe-1907"
#> [19] "TMe-3286" "TMe-1444" "TMe-2860" "TMe-1107" "TMe-3771" "TMe-2866"
#> [25] "TMe-3200" "TMe-2056" "TMe-2952" "TMe-2329" "TMe-3237" "TMe-1385"
#> [31] "TMe-86"  
#> 
#> [[3]]
#>  [1] "TMe-70"   "TMe-1897" "TMe-3383" "TMe-3133" "TMe-2481" "TMe-3445"
#>  [7] "TMe-234"  "TMe-3596" "TMe-1863" "TMe-773"  "TMe-2926" "TMe-116" 
#> [13] "TMe-878"  "TMe-3274" "TMe-200"  "TMe-35"   "TMe-3565" "TMe-1868"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1434" "TMe-3527" "TMe-3255"
#>  [7] "TMe-3542" "TMe-650"  "TMe-2971" "TMe-2809" "TMe-663"  "TMe-1419"
#> [13] "TMe-1511" "TMe-3054" "TMe-3442" "TMe-1988" "TMe-241"  "TMe-2924"
#> [19] "TMe-1123" "TMe-2567" "TMe-3225" "TMe-3730" "TMe-3242" "TMe-550" 
#> [25] "TMe-1020" "TMe-225"  "TMe-2240" "TMe-3422" "TMe-1525" "TMe-3273"
#> [31] "TMe-1376" "TMe-2956" "TMe-3167" "TMe-388"  "TMe-3760"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2907" "TMe-3185" "TMe-277"  "TMe-1381"
#>  [7] "TMe-378"  "TMe-2953" "TMe-723"  "TMe-2589" "TMe-2531" "TMe-1257"
#> [13] "TMe-2853" "TMe-2121" "TMe-2041" "TMe-969"  "TMe-3034" "TMe-712" 
#> [19] "TMe-287"  "TMe-2790" "TMe-1500" "TMe-3356" "TMe-2355" "TMe-284" 
#> [25] "TMe-1004" "TMe-48"   "TMe-2820" "TMe-336"  "TMe-294"  "TMe-3628"
#> [31] "TMe-158"  "TMe-536"  "TMe-1458" "TMe-759"  "TMe-2612" "TMe-1963"
#> [37] "TMe-3311" "TMe-2296" "TMe-208"  "TMe-1357"
#> 
#> [[6]]
#>  [1] "TMe-3549" "TMe-1079" "TMe-1566" "TMe-1900" "TMe-1232" "TMe-1509"
#>  [7] "TMe-693"  "TMe-1883" "TMe-2791" "TMe-3095" "TMe-2035" "TMe-580" 
#> [13] "TMe-1875" "TMe-1428" "TMe-1174" "TMe-751"  "TMe-2064"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Pooled Brillouin diversity index")


# Mean Margalef's richness index
greedysel_best_mean_margalef <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_margalef,
                   metric = "mean", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_mean_margalef
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-3424" "TMe-778"  "TMe-2967"
#>  [7] "TMe-2965" "TMe-841"  "TMe-2010" "TMe-3030" "TMe-1140" "TMe-1915"
#> [13] "TMe-132"  "TMe-1226" "TMe-1360" "TMe-1564" "TMe-3262" "TMe-41"  
#> [19] "TMe-44"   "TMe-117"  "TMe-290"  "TMe-299"  "TMe-300"  "TMe-407" 
#> [25] "TMe-410"  "TMe-438"  "TMe-446"  "TMe-469"  "TMe-489"  "TMe-500" 
#> [31] "TMe-501" 
#> 
#> [[2]]
#>  [1] "TMe-3617" "TMe-3466" "TMe-1619" "TMe-3605" "TMe-1907" "TMe-1385"
#>  [7] "TMe-40"   "TMe-3338" "TMe-2304" "TMe-2033" "TMe-160"  "TMe-289" 
#> [13] "TMe-369"  "TMe-509"  "TMe-1241" "TMe-1461" "TMe-3200" "TMe-3222"
#> [19] "TMe-3237" "TMe-2"    "TMe-39"   "TMe-60"   "TMe-85"   "TMe-339" 
#> [25] "TMe-365"  "TMe-377"  "TMe-404"  "TMe-409"  "TMe-433"  "TMe-450" 
#> [31] "TMe-455" 
#> 
#> [[3]]
#>  [1] "TMe-2481" "TMe-3596" "TMe-200"  "TMe-3383" "TMe-1897" "TMe-3133"
#>  [7] "TMe-405"  "TMe-4"    "TMe-1993" "TMe-123"  "TMe-267"  "TMe-3046"
#> [13] "TMe-1868" "TMe-2592" "TMe-35"   "TMe-2926" "TMe-6"    "TMe-13"  
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1350" "TMe-552"  "TMe-1434"
#>  [7] "TMe-3527" "TMe-3054" "TMe-3255" "TMe-540"  "TMe-3428" "TMe-241" 
#> [13] "TMe-103"  "TMe-266"  "TMe-427"  "TMe-460"  "TMe-550"  "TMe-2809"
#> [19] "TMe-3273" "TMe-1020" "TMe-1511" "TMe-2567" "TMe-11"   "TMe-12"  
#> [25] "TMe-57"   "TMe-62"   "TMe-90"   "TMe-107"  "TMe-108"  "TMe-150" 
#> [31] "TMe-154"  "TMe-170"  "TMe-182"  "TMe-191"  "TMe-204" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2307" "TMe-3185" "TMe-197"  "TMe-723" 
#>  [7] "TMe-378"  "TMe-2953" "TMe-277"  "TMe-2853" "TMe-1004" "TMe-336" 
#> [13] "TMe-399"  "TMe-457"  "TMe-1904" "TMe-2041" "TMe-181"  "TMe-48"  
#> [19] "TMe-98"   "TMe-99"   "TMe-158"  "TMe-167"  "TMe-212"  "TMe-275" 
#> [25] "TMe-294"  "TMe-305"  "TMe-307"  "TMe-332"  "TMe-334"  "TMe-340" 
#> [31] "TMe-347"  "TMe-359"  "TMe-360"  "TMe-362"  "TMe-363"  "TMe-373" 
#> [37] "TMe-376"  "TMe-382"  "TMe-385"  "TMe-390" 
#> 
#> [[6]]
#>  [1] "TMe-1721" "TMe-3095" "TMe-1232" "TMe-2035" "TMe-1428" "TMe-1858"
#>  [7] "TMe-315"  "TMe-693"  "TMe-1875" "TMe-1174" "TMe-130"  "TMe-751" 
#> [13] "TMe-1775" "TMe-856"  "TMe-138"  "TMe-269"  "TMe-321" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Mean Margalef's diversity index")


# Pooled Margalef's richness index
greedysel_best_sum_margalef <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_margalef,
                   metric = "pooled", search = "greedy",
                   local.search = "best.improvement",max.iter = 3)
greedysel_best_sum_margalef
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-677"  "TMe-3705" "TMe-1717"
#>  [7] "TMe-3030" "TMe-1425" "TMe-2010" "TMe-2069" "TMe-1140" "TMe-3396"
#> [13] "TMe-1226" "TMe-778"  "TMe-1360" "TMe-132"  "TMe-1922" "TMe-1915"
#> [19] "TMe-41"   "TMe-44"   "TMe-117"  "TMe-290"  "TMe-299"  "TMe-300" 
#> [25] "TMe-407"  "TMe-410"  "TMe-438"  "TMe-446"  "TMe-469"  "TMe-489" 
#> [31] "TMe-500" 
#> 
#> [[2]]
#>  [1] "TMe-2935" "TMe-2352" "TMe-1107" "TMe-3466" "TMe-1459" "TMe-369" 
#>  [7] "TMe-509"  "TMe-3771" "TMe-47"   "TMe-160"  "TMe-1385" "TMe-1907"
#> [13] "TMe-2033" "TMe-289"  "TMe-3200" "TMe-2304" "TMe-1461" "TMe-3222"
#> [19] "TMe-3237" "TMe-2"    "TMe-39"   "TMe-60"   "TMe-85"   "TMe-339" 
#> [25] "TMe-365"  "TMe-377"  "TMe-404"  "TMe-409"  "TMe-433"  "TMe-450" 
#> [31] "TMe-455" 
#> 
#> [[3]]
#>  [1] "TMe-10"   "TMe-6"    "TMe-3556" "TMe-3715" "TMe-1838" "TMe-3804"
#>  [7] "TMe-234"  "TMe-1897" "TMe-3383" "TMe-35"   "TMe-4"    "TMe-123" 
#> [13] "TMe-267"  "TMe-1993" "TMe-1868" "TMe-2926" "TMe-3596" "TMe-13"  
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1350" "TMe-552"  "TMe-1434"
#>  [7] "TMe-3255" "TMe-3527" "TMe-3054" "TMe-540"  "TMe-3428" "TMe-241" 
#> [13] "TMe-550"  "TMe-103"  "TMe-427"  "TMe-266"  "TMe-460"  "TMe-2809"
#> [19] "TMe-3273" "TMe-1020" "TMe-1511" "TMe-2567" "TMe-11"   "TMe-12"  
#> [25] "TMe-57"   "TMe-62"   "TMe-90"   "TMe-107"  "TMe-108"  "TMe-150" 
#> [31] "TMe-154"  "TMe-170"  "TMe-182"  "TMe-191"  "TMe-204" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2307" "TMe-3185" "TMe-197"  "TMe-723" 
#>  [7] "TMe-378"  "TMe-2953" "TMe-277"  "TMe-2853" "TMe-1004" "TMe-336" 
#> [13] "TMe-399"  "TMe-457"  "TMe-1904" "TMe-2041" "TMe-181"  "TMe-48"  
#> [19] "TMe-98"   "TMe-99"   "TMe-158"  "TMe-167"  "TMe-212"  "TMe-275" 
#> [25] "TMe-294"  "TMe-305"  "TMe-307"  "TMe-332"  "TMe-334"  "TMe-340" 
#> [31] "TMe-347"  "TMe-359"  "TMe-360"  "TMe-362"  "TMe-363"  "TMe-373" 
#> [37] "TMe-376"  "TMe-382"  "TMe-385"  "TMe-390" 
#> 
#> [[6]]
#>  [1] "TMe-1592" "TMe-1416" "TMe-693"  "TMe-1554" "TMe-315"  "TMe-2067"
#>  [7] "TMe-1883" "TMe-1124" "TMe-751"  "TMe-1875" "TMe-130"  "TMe-1858"
#> [13] "TMe-856"  "TMe-1775" "TMe-963"  "TMe-3095" "TMe-138" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt best improvement",
       subtitle = "Pooled Margalef's diversity index")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Greedy search with 1-opt first improvement
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Mean richness
greedysel_first_mean_richness <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "richness",
                   metric = "mean", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_mean_richness
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-677"  "TMe-3475" "TMe-1717"
#>  [7] "TMe-841"  "TMe-1425" "TMe-2010" "TMe-3030" "TMe-3292" "TMe-2069"
#> [13] "TMe-132"  "TMe-778"  "TMe-1226" "TMe-1360" "TMe-1922" "TMe-3262"
#> [19] "TMe-1915" "TMe-41"   "TMe-44"   "TMe-117"  "TMe-290"  "TMe-299" 
#> [25] "TMe-300"  "TMe-407"  "TMe-410"  "TMe-438"  "TMe-446"  "TMe-469" 
#> [31] "TMe-489" 
#> 
#> [[2]]
#>  [1] "TMe-160"  "TMe-1107" "TMe-2127" "TMe-1444" "TMe-1698" "TMe-1907"
#>  [7] "TMe-3466" "TMe-3338" "TMe-289"  "TMe-2033" "TMe-3200" "TMe-369" 
#> [13] "TMe-509"  "TMe-1241" "TMe-1385" "TMe-2304" "TMe-1461" "TMe-3222"
#> [19] "TMe-3237" "TMe-2"    "TMe-39"   "TMe-60"   "TMe-85"   "TMe-339" 
#> [25] "TMe-365"  "TMe-377"  "TMe-404"  "TMe-409"  "TMe-433"  "TMe-450" 
#> [31] "TMe-455" 
#> 
#> [[3]]
#>  [1] "TMe-3216" "TMe-2926" "TMe-785"  "TMe-946"  "TMe-2977" "TMe-116" 
#>  [7] "TMe-420"  "TMe-773"  "TMe-1993" "TMe-3596" "TMe-123"  "TMe-267" 
#> [13] "TMe-1868" "TMe-2592" "TMe-3046" "TMe-35"   "TMe-4"    "TMe-6"   
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1350" "TMe-552"  "TMe-1434"
#>  [7] "TMe-3255" "TMe-3527" "TMe-3054" "TMe-540"  "TMe-3428" "TMe-241" 
#> [13] "TMe-103"  "TMe-266"  "TMe-427"  "TMe-460"  "TMe-550"  "TMe-2809"
#> [19] "TMe-3273" "TMe-1020" "TMe-1511" "TMe-2567" "TMe-11"   "TMe-12"  
#> [25] "TMe-57"   "TMe-62"   "TMe-90"   "TMe-107"  "TMe-108"  "TMe-150" 
#> [31] "TMe-154"  "TMe-170"  "TMe-182"  "TMe-191"  "TMe-204" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2307" "TMe-3185" "TMe-197"  "TMe-723" 
#>  [7] "TMe-378"  "TMe-2953" "TMe-277"  "TMe-2853" "TMe-1004" "TMe-336" 
#> [13] "TMe-399"  "TMe-457"  "TMe-1904" "TMe-2041" "TMe-181"  "TMe-48"  
#> [19] "TMe-98"   "TMe-99"   "TMe-158"  "TMe-167"  "TMe-212"  "TMe-275" 
#> [25] "TMe-294"  "TMe-305"  "TMe-307"  "TMe-332"  "TMe-334"  "TMe-340" 
#> [31] "TMe-347"  "TMe-359"  "TMe-360"  "TMe-362"  "TMe-363"  "TMe-373" 
#> [37] "TMe-376"  "TMe-382"  "TMe-385"  "TMe-390" 
#> 
#> [[6]]
#>  [1] "TMe-268"  "TMe-1053" "TMe-1566" "TMe-1174" "TMe-2791" "TMe-2551"
#>  [7] "TMe-1883" "TMe-315"  "TMe-751"  "TMe-693"  "TMe-2510" "TMe-465" 
#> [13] "TMe-856"  "TMe-3095" "TMe-130"  "TMe-138"  "TMe-269" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Mean richness")


# Pooled richness
greedysel_first_sum_richness <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "richness",
                   metric = "pooled", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_sum_richness
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-677"  "TMe-3475" "TMe-1717"
#>  [7] "TMe-841"  "TMe-1425" "TMe-2010" "TMe-3030" "TMe-3292" "TMe-2069"
#> [13] "TMe-132"  "TMe-778"  "TMe-1226" "TMe-1360" "TMe-1922" "TMe-3262"
#> [19] "TMe-1915" "TMe-41"   "TMe-44"   "TMe-117"  "TMe-290"  "TMe-299" 
#> [25] "TMe-300"  "TMe-407"  "TMe-410"  "TMe-438"  "TMe-446"  "TMe-469" 
#> [31] "TMe-489" 
#> 
#> [[2]]
#>  [1] "TMe-3690" "TMe-509"  "TMe-1919" "TMe-2127" "TMe-2851" "TMe-3141"
#>  [7] "TMe-3466" "TMe-3338" "TMe-160"  "TMe-2033" "TMe-2860" "TMe-369" 
#> [13] "TMe-1241" "TMe-2304" "TMe-1461" "TMe-3200" "TMe-3222" "TMe-3237"
#> [19] "TMe-2"    "TMe-39"   "TMe-60"   "TMe-85"   "TMe-289"  "TMe-339" 
#> [25] "TMe-365"  "TMe-377"  "TMe-404"  "TMe-409"  "TMe-433"  "TMe-450" 
#> [31] "TMe-455" 
#> 
#> [[3]]
#>  [1] "TMe-3592" "TMe-1176" "TMe-3274" "TMe-946"  "TMe-1910" "TMe-3565"
#>  [7] "TMe-35"   "TMe-1593" "TMe-2926" "TMe-381"  "TMe-123"  "TMe-267" 
#> [13] "TMe-1868" "TMe-1993" "TMe-2592" "TMe-3046" "TMe-3596" "TMe-4"   
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1350" "TMe-552"  "TMe-1434"
#>  [7] "TMe-3255" "TMe-3527" "TMe-3054" "TMe-540"  "TMe-3428" "TMe-241" 
#> [13] "TMe-103"  "TMe-266"  "TMe-427"  "TMe-460"  "TMe-550"  "TMe-2809"
#> [19] "TMe-3273" "TMe-1020" "TMe-1511" "TMe-2567" "TMe-11"   "TMe-12"  
#> [25] "TMe-57"   "TMe-62"   "TMe-90"   "TMe-107"  "TMe-108"  "TMe-150" 
#> [31] "TMe-154"  "TMe-170"  "TMe-182"  "TMe-191"  "TMe-204" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2307" "TMe-3185" "TMe-197"  "TMe-723" 
#>  [7] "TMe-378"  "TMe-2953" "TMe-277"  "TMe-2853" "TMe-1004" "TMe-336" 
#> [13] "TMe-399"  "TMe-457"  "TMe-1904" "TMe-2041" "TMe-181"  "TMe-48"  
#> [19] "TMe-98"   "TMe-99"   "TMe-158"  "TMe-167"  "TMe-212"  "TMe-275" 
#> [25] "TMe-294"  "TMe-305"  "TMe-307"  "TMe-332"  "TMe-334"  "TMe-340" 
#> [31] "TMe-347"  "TMe-359"  "TMe-360"  "TMe-362"  "TMe-363"  "TMe-373" 
#> [37] "TMe-376"  "TMe-382"  "TMe-385"  "TMe-390" 
#> 
#> [[6]]
#>  [1] "TMe-1007" "TMe-3095" "TMe-1232" "TMe-1428" "TMe-1858" "TMe-2035"
#>  [7] "TMe-315"  "TMe-751"  "TMe-1875" "TMe-693"  "TMe-1174" "TMe-1775"
#> [13] "TMe-130"  "TMe-856"  "TMe-138"  "TMe-269"  "TMe-321" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Pooled richness")


# Mean Shannon-Weaver diversity index
greedysel_first_mean_shannon <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "shannon",
                   metric = "mean", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_mean_shannon
#> [[1]]
#>  [1] "TMe-3423" "TMe-3394" "TMe-867"  "TMe-2027" "TMe-2069" "TMe-1425"
#>  [7] "TMe-41"   "TMe-2010" "TMe-3292" "TMe-566"  "TMe-3030" "TMe-1360"
#> [13] "TMe-2993" "TMe-2965" "TMe-1960" "TMe-2967" "TMe-407"  "TMe-2004"
#> [19] "TMe-606"  "TMe-1190" "TMe-3534" "TMe-3396" "TMe-2785" "TMe-2237"
#> [25] "TMe-3112" "TMe-1922" "TMe-3102" "TMe-1096" "TMe-778"  "TMe-1786"
#> [31] "TMe-3475"
#> 
#> [[2]]
#>  [1] "TMe-86"   "TMe-1444" "TMe-2056" "TMe-47"   "TMe-796"  "TMe-2952"
#>  [7] "TMe-2033" "TMe-2860" "TMe-3771" "TMe-478"  "TMe-3466" "TMe-3338"
#> [13] "TMe-3605" "TMe-2866" "TMe-509"  "TMe-3141" "TMe-3237" "TMe-369" 
#> [19] "TMe-1919" "TMe-3200" "TMe-40"   "TMe-1907" "TMe-3009" "TMe-1241"
#> [25] "TMe-3766" "TMe-1461" "TMe-160"  "TMe-2128" "TMe-1107" "TMe-3286"
#> [31] "TMe-289" 
#> 
#> [[3]]
#>  [1] "TMe-200"  "TMe-3596" "TMe-1200" "TMe-3565" "TMe-3274" "TMe-946" 
#>  [7] "TMe-2823" "TMe-1868" "TMe-1993" "TMe-35"   "TMe-405"  "TMe-2926"
#> [13] "TMe-3383" "TMe-381"  "TMe-1897" "TMe-2481" "TMe-234"  "TMe-3715"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-3225" "TMe-1434" "TMe-3255"
#>  [7] "TMe-3267" "TMe-1020" "TMe-3730" "TMe-3167" "TMe-225"  "TMe-2809"
#> [13] "TMe-1988" "TMe-1511" "TMe-3527" "TMe-3442" "TMe-3054" "TMe-1330"
#> [19] "TMe-2971" "TMe-241"  "TMe-2924" "TMe-1376" "TMe-2956" "TMe-3273"
#> [25] "TMe-2788" "TMe-266"  "TMe-1525" "TMe-3242" "TMe-3542" "TMe-550" 
#> [31] "TMe-2567" "TMe-3760" "TMe-1123" "TMe-1419" "TMe-663" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1534" "TMe-3185" "TMe-197"  "TMe-3034"
#>  [7] "TMe-277"  "TMe-2953" "TMe-378"  "TMe-723"  "TMe-2041" "TMe-2853"
#> [13] "TMe-2589" "TMe-2531" "TMe-284"  "TMe-1381" "TMe-1458" "TMe-969" 
#> [19] "TMe-712"  "TMe-1004" "TMe-2121" "TMe-2355" "TMe-2790" "TMe-287" 
#> [25] "TMe-336"  "TMe-3356" "TMe-1500" "TMe-536"  "TMe-2820" "TMe-3332"
#> [31] "TMe-1963" "TMe-294"  "TMe-3628" "TMe-759"  "TMe-2612" "TMe-48"  
#> [37] "TMe-1257" "TMe-2863" "TMe-1099" "TMe-2907"
#> 
#> [[6]]
#>  [1] "TMe-1775" "TMe-2983" "TMe-2064" "TMe-3095" "TMe-856"  "TMe-1174"
#>  [7] "TMe-1883" "TMe-1232" "TMe-693"  "TMe-1509" "TMe-1858" "TMe-751" 
#> [13] "TMe-3549" "TMe-2067" "TMe-1875" "TMe-315"  "TMe-1428"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Mean Shannon-Weaver diversity index")


# Pooled Shannon-Weaver diversity index
greedysel_first_sum_shannon <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "shannon",
                   metric = "pooled", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_sum_shannon
#> [[1]]
#>  [1] "TMe-3423" "TMe-3394" "TMe-867"  "TMe-2027" "TMe-2069" "TMe-1425"
#>  [7] "TMe-41"   "TMe-2010" "TMe-3292" "TMe-566"  "TMe-3030" "TMe-1360"
#> [13] "TMe-2993" "TMe-2965" "TMe-1960" "TMe-2967" "TMe-407"  "TMe-2004"
#> [19] "TMe-606"  "TMe-1190" "TMe-3534" "TMe-3396" "TMe-2785" "TMe-2237"
#> [25] "TMe-3112" "TMe-1922" "TMe-3102" "TMe-1096" "TMe-778"  "TMe-1786"
#> [31] "TMe-3475"
#> 
#> [[2]]
#>  [1] "TMe-3605" "TMe-2"    "TMe-1107" "TMe-1444" "TMe-3237" "TMe-3466"
#>  [7] "TMe-1907" "TMe-509"  "TMe-47"   "TMe-3141" "TMe-1919" "TMe-289" 
#> [13] "TMe-369"  "TMe-2866" "TMe-3200" "TMe-2033" "TMe-3338" "TMe-3771"
#> [19] "TMe-2056" "TMe-2952" "TMe-1461" "TMe-3805" "TMe-2128" "TMe-160" 
#> [25] "TMe-478"  "TMe-2860" "TMe-3286" "TMe-40"   "TMe-3009" "TMe-3690"
#> [31] "TMe-1241"
#> 
#> [[3]]
#>  [1] "TMe-3596" "TMe-200"  "TMe-1838" "TMe-1593" "TMe-1897" "TMe-267" 
#>  [7] "TMe-3128" "TMe-141"  "TMe-3715" "TMe-785"  "TMe-2926" "TMe-2481"
#> [13] "TMe-3383" "TMe-116"  "TMe-3565" "TMe-1868" "TMe-35"   "TMe-3274"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-3225" "TMe-1434" "TMe-3255"
#>  [7] "TMe-3267" "TMe-1020" "TMe-3730" "TMe-3167" "TMe-225"  "TMe-2809"
#> [13] "TMe-1988" "TMe-1511" "TMe-3527" "TMe-3442" "TMe-3054" "TMe-1330"
#> [19] "TMe-2971" "TMe-241"  "TMe-2924" "TMe-1376" "TMe-2956" "TMe-3273"
#> [25] "TMe-2788" "TMe-266"  "TMe-1525" "TMe-3242" "TMe-3542" "TMe-550" 
#> [31] "TMe-2567" "TMe-3760" "TMe-1123" "TMe-1419" "TMe-663" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1534" "TMe-3185" "TMe-197"  "TMe-3034"
#>  [7] "TMe-277"  "TMe-2953" "TMe-378"  "TMe-723"  "TMe-2041" "TMe-2853"
#> [13] "TMe-2589" "TMe-2531" "TMe-284"  "TMe-1381" "TMe-1458" "TMe-969" 
#> [19] "TMe-712"  "TMe-1004" "TMe-2121" "TMe-2355" "TMe-2790" "TMe-287" 
#> [25] "TMe-336"  "TMe-3356" "TMe-1500" "TMe-536"  "TMe-2820" "TMe-3332"
#> [31] "TMe-1963" "TMe-294"  "TMe-3628" "TMe-759"  "TMe-2612" "TMe-48"  
#> [37] "TMe-1257" "TMe-2863" "TMe-1099" "TMe-2907"
#> 
#> [[6]]
#>  [1] "TMe-1353" "TMe-3095" "TMe-2791" "TMe-1756" "TMe-1124" "TMe-1509"
#>  [7] "TMe-1539" "TMe-751"  "TMe-1875" "TMe-2510" "TMe-693"  "TMe-3549"
#> [13] "TMe-580"  "TMe-2064" "TMe-1566" "TMe-1174" "TMe-130" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Pooled Shannon-Weaver diversity index")


# Mean Gini-Simpson diversity index
greedysel_first_mean_simpson <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "simpson",
                   metric = "mean", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_mean_simpson
#> [[1]]
#>  [1] "TMe-3423" "TMe-2785" "TMe-867"  "TMe-3452" "TMe-2069" "TMe-1425"
#>  [7] "TMe-3389" "TMe-1890" "TMe-41"   "TMe-3292" "TMe-1360" "TMe-1960"
#> [13] "TMe-2965" "TMe-1190" "TMe-2967" "TMe-3030" "TMe-1423" "TMe-28"  
#> [19] "TMe-606"  "TMe-3475" "TMe-2993" "TMe-3102" "TMe-2372" "TMe-841" 
#> [25] "TMe-3319" "TMe-132"  "TMe-3208" "TMe-1086" "TMe-3705" "TMe-1922"
#> [31] "TMe-1564"
#> 
#> [[2]]
#>  [1] "TMe-636"  "TMe-160"  "TMe-1107" "TMe-3141" "TMe-3338" "TMe-2127"
#>  [7] "TMe-3466" "TMe-509"  "TMe-674"  "TMe-47"   "TMe-1385" "TMe-1754"
#> [13] "TMe-86"   "TMe-369"  "TMe-1919" "TMe-3690" "TMe-2033" "TMe-2056"
#> [19] "TMe-1444" "TMe-2866" "TMe-3200" "TMe-3805" "TMe-1907" "TMe-3771"
#> [25] "TMe-2952" "TMe-171"  "TMe-2128" "TMe-3140" "TMe-40"   "TMe-2329"
#> [31] "TMe-2860"
#> 
#> [[3]]
#>  [1] "TMe-141"  "TMe-785"  "TMe-2756" "TMe-2926" "TMe-773"  "TMe-3274"
#>  [7] "TMe-3565" "TMe-1939" "TMe-3596" "TMe-1897" "TMe-3383" "TMe-2481"
#> [13] "TMe-200"  "TMe-267"  "TMe-3143" "TMe-1200" "TMe-2897" "TMe-1868"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1434" "TMe-388"  "TMe-3275"
#>  [7] "TMe-1083" "TMe-3225" "TMe-2809" "TMe-2924" "TMe-1020" "TMe-962" 
#> [13] "TMe-812"  "TMe-1511" "TMe-3442" "TMe-3054" "TMe-2971" "TMe-1988"
#> [19] "TMe-241"  "TMe-2552" "TMe-3417" "TMe-3542" "TMe-3527" "TMe-1525"
#> [25] "TMe-650"  "TMe-180"  "TMe-2240" "TMe-1419" "TMe-550"  "TMe-3242"
#> [31] "TMe-3760" "TMe-3198" "TMe-2956" "TMe-3273" "TMe-3390"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1399" "TMe-3185" "TMe-712"  "TMe-1381"
#>  [7] "TMe-2589" "TMe-3736" "TMe-723"  "TMe-1458" "TMe-2953" "TMe-2531"
#> [13] "TMe-378"  "TMe-284"  "TMe-1534" "TMe-3034" "TMe-294"  "TMe-969" 
#> [19] "TMe-2355" "TMe-1500" "TMe-3356" "TMe-2820" "TMe-536"  "TMe-2041"
#> [25] "TMe-759"  "TMe-2003" "TMe-2853" "TMe-2121" "TMe-3311" "TMe-1099"
#> [31] "TMe-1357" "TMe-2612" "TMe-208"  "TMe-336"  "TMe-158"  "TMe-2907"
#> [37] "TMe-2790" "TMe-3628" "TMe-48"   "TMe-277" 
#> 
#> [[6]]
#>  [1] "TMe-130"  "TMe-3095" "TMe-2983" "TMe-1232" "TMe-1353" "TMe-1035"
#>  [7] "TMe-2791" "TMe-3549" "TMe-1875" "TMe-693"  "TMe-1509" "TMe-465" 
#> [13] "TMe-3116" "TMe-1539" "TMe-580"  "TMe-725"  "TMe-1428"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Mean Gini-Simpson diversity index")


# Pooled Gini-Simpson diversity index
greedysel_first_sum_simpson <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "simpson",
                   metric = "pooled", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_sum_simpson
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-3282" "TMe-41"   "TMe-839"  "TMe-1960"
#>  [7] "TMe-2965" "TMe-2967" "TMe-1423" "TMe-1190" "TMe-3292" "TMe-3437"
#> [13] "TMe-1425" "TMe-2785" "TMe-867"  "TMe-3394" "TMe-3102" "TMe-2993"
#> [19] "TMe-1360" "TMe-3475" "TMe-841"  "TMe-2462" "TMe-28"   "TMe-2069"
#> [25] "TMe-606"  "TMe-717"  "TMe-1600" "TMe-3151" "TMe-3319" "TMe-1786"
#> [31] "TMe-2927"
#> 
#> [[2]]
#>  [1] "TMe-3009" "TMe-47"   "TMe-2056" "TMe-86"   "TMe-1444" "TMe-1107"
#>  [7] "TMe-3141" "TMe-3338" "TMe-2866" "TMe-3466" "TMe-369"  "TMe-2033"
#> [13] "TMe-1754" "TMe-3200" "TMe-509"  "TMe-1907" "TMe-1919" "TMe-3237"
#> [19] "TMe-3771" "TMe-2952" "TMe-3766" "TMe-2128" "TMe-1385" "TMe-40"  
#> [25] "TMe-2860" "TMe-2329" "TMe-3805" "TMe-478"  "TMe-160"  "TMe-3286"
#> [31] "TMe-3140"
#> 
#> [[3]]
#>  [1] "TMe-773"  "TMe-3445" "TMe-2926" "TMe-878"  "TMe-116"  "TMe-3274"
#>  [7] "TMe-1897" "TMe-3383" "TMe-3596" "TMe-200"  "TMe-2481" "TMe-141" 
#> [13] "TMe-785"  "TMe-3715" "TMe-1593" "TMe-1868" "TMe-267"  "TMe-3565"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1434" "TMe-388"  "TMe-3275"
#>  [7] "TMe-1406" "TMe-3225" "TMe-2809" "TMe-2924" "TMe-1020" "TMe-1511"
#> [13] "TMe-3707" "TMe-2971" "TMe-3442" "TMe-550"  "TMe-225"  "TMe-241" 
#> [19] "TMe-2956" "TMe-2552" "TMe-824"  "TMe-3730" "TMe-3527" "TMe-1525"
#> [25] "TMe-3542" "TMe-3054" "TMe-1328" "TMe-1988" "TMe-3167" "TMe-3417"
#> [31] "TMe-180"  "TMe-3273" "TMe-3390" "TMe-3576" "TMe-3242"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-1399" "TMe-3185" "TMe-712"  "TMe-1381"
#>  [7] "TMe-2589" "TMe-3736" "TMe-723"  "TMe-1458" "TMe-2953" "TMe-2531"
#> [13] "TMe-378"  "TMe-284"  "TMe-1534" "TMe-3034" "TMe-1357" "TMe-2820"
#> [19] "TMe-1500" "TMe-3356" "TMe-294"  "TMe-3311" "TMe-2853" "TMe-2121"
#> [25] "TMe-2041" "TMe-158"  "TMe-536"  "TMe-969"  "TMe-759"  "TMe-208" 
#> [31] "TMe-2003" "TMe-1099" "TMe-2907" "TMe-336"  "TMe-3628" "TMe-2612"
#> [37] "TMe-48"   "TMe-2790" "TMe-2355" "TMe-277" 
#> 
#> [[6]]
#>  [1] "TMe-2064" "TMe-836"  "TMe-1769" "TMe-1633" "TMe-693"  "TMe-3095"
#>  [7] "TMe-361"  "TMe-1035" "TMe-1232" "TMe-1428" "TMe-2983" "TMe-1875"
#> [13] "TMe-3549" "TMe-2791" "TMe-1678" "TMe-1353" "TMe-1509"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Pooled Gini-Simpson diversity index")


# Mean McIntosh diversity index
greedysel_first_mean_mcintosh <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "mcintosh",
                   metric = "pooled", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_mean_mcintosh
#> [[1]]
#>  [1] "TMe-3423" "TMe-3208" "TMe-867"  "TMe-132"  "TMe-3359" "TMe-2967"
#>  [7] "TMe-1425" "TMe-2785" "TMe-1190" "TMe-3437" "TMe-670"  "TMe-2934"
#> [13] "TMe-2993" "TMe-1360" "TMe-2965" "TMe-3433" "TMe-407"  "TMe-3292"
#> [19] "TMe-41"   "TMe-1960" "TMe-2069" "TMe-1423" "TMe-3102" "TMe-606" 
#> [25] "TMe-3396" "TMe-28"   "TMe-3236" "TMe-1333" "TMe-3319" "TMe-3030"
#> [31] "TMe-3705"
#> 
#> [[2]]
#>  [1] "TMe-1241" "TMe-3338" "TMe-160"  "TMe-1107" "TMe-2"    "TMe-3466"
#>  [7] "TMe-509"  "TMe-3141" "TMe-1919" "TMe-369"  "TMe-1907" "TMe-1754"
#> [13] "TMe-2128" "TMe-3200" "TMe-47"   "TMe-1385" "TMe-2056" "TMe-2952"
#> [19] "TMe-3766" "TMe-2033" "TMe-3771" "TMe-2866" "TMe-3690" "TMe-1444"
#> [25] "TMe-2860" "TMe-3286" "TMe-3237" "TMe-2329" "TMe-636"  "TMe-3009"
#> [31] "TMe-86"  
#> 
#> [[3]]
#>  [1] "TMe-773"  "TMe-1939" "TMe-3274" "TMe-1897" "TMe-946"  "TMe-3596"
#>  [7] "TMe-785"  "TMe-3383" "TMe-13"   "TMe-2926" "TMe-261"  "TMe-1593"
#> [13] "TMe-405"  "TMe-1993" "TMe-141"  "TMe-3445" "TMe-2481" "TMe-2823"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-2552" "TMe-1434" "TMe-3255"
#>  [7] "TMe-3442" "TMe-2971" "TMe-209"  "TMe-2809" "TMe-1575" "TMe-1377"
#> [13] "TMe-2924" "TMe-3707" "TMe-1511" "TMe-241"  "TMe-388"  "TMe-2956"
#> [19] "TMe-550"  "TMe-3225" "TMe-650"  "TMe-3422" "TMe-2755" "TMe-225" 
#> [25] "TMe-1020" "TMe-1597" "TMe-3054" "TMe-3527" "TMe-3542" "TMe-1525"
#> [31] "TMe-1376" "TMe-3417" "TMe-3390" "TMe-1988" "TMe-2240"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2612" "TMe-3185" "TMe-2790" "TMe-1381"
#>  [7] "TMe-378"  "TMe-2953" "TMe-723"  "TMe-2589" "TMe-2531" "TMe-344" 
#> [13] "TMe-2121" "TMe-2853" "TMe-712"  "TMe-3332" "TMe-1924" "TMe-3034"
#> [19] "TMe-1500" "TMe-1357" "TMe-208"  "TMe-969"  "TMe-287"  "TMe-284" 
#> [25] "TMe-3356" "TMe-294"  "TMe-759"  "TMe-536"  "TMe-2355" "TMe-2820"
#> [31] "TMe-158"  "TMe-2041" "TMe-2750" "TMe-3628" "TMe-1458" "TMe-3311"
#> [37] "TMe-2003" "TMe-1099" "TMe-2863" "TMe-755" 
#> 
#> [[6]]
#>  [1] "TMe-361"  "TMe-2510" "TMe-1353" "TMe-3549" "TMe-1875" "TMe-693" 
#>  [7] "TMe-580"  "TMe-2064" "TMe-3095" "TMe-1035" "TMe-751"  "TMe-1428"
#> [13] "TMe-1232" "TMe-2035" "TMe-2791" "TMe-236"  "TMe-470" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Mean McIntosh diversity index")


# Pooled McIntosh diversity index
greedysel_first_sum_mcintosh <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.index = "mcintosh",
                   metric = "pooled", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_sum_mcintosh
#> [[1]]
#>  [1] "TMe-3423" "TMe-3208" "TMe-867"  "TMe-132"  "TMe-3359" "TMe-2967"
#>  [7] "TMe-1425" "TMe-2785" "TMe-1190" "TMe-3437" "TMe-670"  "TMe-2934"
#> [13] "TMe-2993" "TMe-1360" "TMe-2965" "TMe-3433" "TMe-407"  "TMe-3292"
#> [19] "TMe-41"   "TMe-1960" "TMe-2069" "TMe-1423" "TMe-3102" "TMe-606" 
#> [25] "TMe-3396" "TMe-28"   "TMe-3236" "TMe-1333" "TMe-3319" "TMe-3030"
#> [31] "TMe-3705"
#> 
#> [[2]]
#>  [1] "TMe-636"  "TMe-2329" "TMe-1385" "TMe-3466" "TMe-369"  "TMe-2056"
#>  [7] "TMe-3338" "TMe-3141" "TMe-47"   "TMe-2866" "TMe-3200" "TMe-509" 
#> [13] "TMe-1907" "TMe-1919" "TMe-1754" "TMe-2033" "TMe-2952" "TMe-3771"
#> [19] "TMe-3766" "TMe-2128" "TMe-3690" "TMe-40"   "TMe-160"  "TMe-1107"
#> [25] "TMe-3805" "TMe-1444" "TMe-2860" "TMe-3286" "TMe-3237" "TMe-478" 
#> [31] "TMe-3140"
#> 
#> [[3]]
#>  [1] "TMe-1868" "TMe-3094" "TMe-3715" "TMe-141"  "TMe-1838" "TMe-3274"
#>  [7] "TMe-773"  "TMe-1897" "TMe-3596" "TMe-3383" "TMe-785"  "TMe-200" 
#> [13] "TMe-2926" "TMe-1593" "TMe-3565" "TMe-878"  "TMe-1993" "TMe-2823"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-2552" "TMe-1434" "TMe-3255"
#>  [7] "TMe-3442" "TMe-2971" "TMe-209"  "TMe-2809" "TMe-1575" "TMe-1377"
#> [13] "TMe-2924" "TMe-3707" "TMe-1511" "TMe-241"  "TMe-388"  "TMe-2956"
#> [19] "TMe-550"  "TMe-3225" "TMe-650"  "TMe-3422" "TMe-2755" "TMe-225" 
#> [25] "TMe-1020" "TMe-1597" "TMe-3054" "TMe-3527" "TMe-3542" "TMe-1525"
#> [31] "TMe-1376" "TMe-3417" "TMe-3390" "TMe-1988" "TMe-2240"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2612" "TMe-3185" "TMe-2790" "TMe-1381"
#>  [7] "TMe-378"  "TMe-2953" "TMe-723"  "TMe-2589" "TMe-2531" "TMe-344" 
#> [13] "TMe-2121" "TMe-2853" "TMe-712"  "TMe-3332" "TMe-1924" "TMe-3034"
#> [19] "TMe-1500" "TMe-1357" "TMe-208"  "TMe-969"  "TMe-287"  "TMe-284" 
#> [25] "TMe-3356" "TMe-294"  "TMe-759"  "TMe-536"  "TMe-2355" "TMe-2820"
#> [31] "TMe-158"  "TMe-2041" "TMe-2750" "TMe-3628" "TMe-1458" "TMe-3311"
#> [37] "TMe-2003" "TMe-1099" "TMe-2863" "TMe-755" 
#> 
#> [[6]]
#>  [1] "TMe-269"  "TMe-1428" "TMe-1232" "TMe-1079" "TMe-3095" "TMe-693" 
#>  [7] "TMe-213"  "TMe-1353" "TMe-3549" "TMe-2791" "TMe-1875" "TMe-470" 
#> [13] "TMe-2067" "TMe-2510" "TMe-2983" "TMe-1403" "TMe-315" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Pooled McIntosh diversity index")


# Mean Brillouin diversity index
greedysel_first_mean_brillouin <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_brillouin,
                   metric = "mean", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_mean_brillouin
#> [[1]]
#>  [1] "TMe-3423" "TMe-3236" "TMe-867"  "TMe-2996" "TMe-2069" "TMe-1425"
#>  [7] "TMe-3030" "TMe-1360" "TMe-3389" "TMe-3292" "TMe-41"   "TMe-2010"
#> [13] "TMe-1190" "TMe-2965" "TMe-407"  "TMe-2785" "TMe-2993" "TMe-1423"
#> [19] "TMe-2967" "TMe-606"  "TMe-2237" "TMe-1960" "TMe-3396" "TMe-1922"
#> [25] "TMe-3102" "TMe-778"  "TMe-132"  "TMe-2934" "TMe-1786" "TMe-3394"
#> [31] "TMe-566" 
#> 
#> [[2]]
#>  [1] "TMe-3009" "TMe-925"  "TMe-47"   "TMe-1385" "TMe-1907" "TMe-477" 
#>  [7] "TMe-509"  "TMe-3466" "TMe-40"   "TMe-3141" "TMe-1241" "TMe-3338"
#> [13] "TMe-160"  "TMe-2128" "TMe-369"  "TMe-3200" "TMe-1919" "TMe-2033"
#> [19] "TMe-1754" "TMe-2860" "TMe-3286" "TMe-3771" "TMe-1461" "TMe-1107"
#> [25] "TMe-3766" "TMe-1444" "TMe-2056" "TMe-2952" "TMe-2866" "TMe-3690"
#> [31] "TMe-3237"
#> 
#> [[3]]
#>  [1] "TMe-2977" "TMe-267"  "TMe-70"   "TMe-3274" "TMe-1804" "TMe-2926"
#>  [7] "TMe-3383" "TMe-1897" "TMe-261"  "TMe-2802" "TMe-1993" "TMe-2823"
#> [13] "TMe-405"  "TMe-2532" "TMe-3596" "TMe-2481" "TMe-234"  "TMe-3445"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1434" "TMe-3527" "TMe-3255"
#>  [7] "TMe-3542" "TMe-650"  "TMe-2971" "TMe-2809" "TMe-663"  "TMe-1419"
#> [13] "TMe-1511" "TMe-3054" "TMe-3442" "TMe-1988" "TMe-241"  "TMe-2924"
#> [19] "TMe-2240" "TMe-2567" "TMe-3225" "TMe-3730" "TMe-3242" "TMe-550" 
#> [25] "TMe-1020" "TMe-225"  "TMe-1123" "TMe-3422" "TMe-1525" "TMe-3273"
#> [31] "TMe-1376" "TMe-2956" "TMe-3167" "TMe-388"  "TMe-3760"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2907" "TMe-3185" "TMe-277"  "TMe-1381"
#>  [7] "TMe-378"  "TMe-2953" "TMe-723"  "TMe-2589" "TMe-2531" "TMe-344" 
#> [13] "TMe-2853" "TMe-2121" "TMe-2041" "TMe-969"  "TMe-3034" "TMe-712" 
#> [19] "TMe-287"  "TMe-2790" "TMe-1500" "TMe-3356" "TMe-2355" "TMe-284" 
#> [25] "TMe-1004" "TMe-48"   "TMe-2820" "TMe-336"  "TMe-294"  "TMe-3628"
#> [31] "TMe-158"  "TMe-536"  "TMe-1458" "TMe-759"  "TMe-1924" "TMe-1963"
#> [37] "TMe-3311" "TMe-2296" "TMe-208"  "TMe-1357"
#> 
#> [[6]]
#>  [1] "TMe-1816" "TMe-465"  "TMe-1053" "TMe-1503" "TMe-3095" "TMe-1232"
#>  [7] "TMe-693"  "TMe-1875" "TMe-1035" "TMe-1428" "TMe-751"  "TMe-580" 
#> [13] "TMe-2064" "TMe-2510" "TMe-1353" "TMe-1174" "TMe-3549"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Mean Brillouin diversity index")


# Pooled Brillouin diversity index
greedysel_first_sum_brillouin <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_brillouin,
                   metric = "pooled", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_sum_brillouin
#> [[1]]
#>  [1] "TMe-3423" "TMe-3236" "TMe-867"  "TMe-2996" "TMe-2069" "TMe-1425"
#>  [7] "TMe-3030" "TMe-1360" "TMe-3389" "TMe-3292" "TMe-41"   "TMe-2010"
#> [13] "TMe-1190" "TMe-2965" "TMe-407"  "TMe-2785" "TMe-2993" "TMe-1423"
#> [19] "TMe-2967" "TMe-606"  "TMe-2237" "TMe-1960" "TMe-3396" "TMe-1922"
#> [25] "TMe-3102" "TMe-778"  "TMe-132"  "TMe-2934" "TMe-1786" "TMe-3394"
#> [31] "TMe-566" 
#> 
#> [[2]]
#>  [1] "TMe-1623" "TMe-2056" "TMe-47"   "TMe-1444" "TMe-3466" "TMe-369" 
#>  [7] "TMe-3338" "TMe-3141" "TMe-1107" "TMe-1907" "TMe-1385" "TMe-1754"
#> [13] "TMe-160"  "TMe-509"  "TMe-3690" "TMe-1919" "TMe-3286" "TMe-2033"
#> [19] "TMe-2860" "TMe-40"   "TMe-3605" "TMe-3771" "TMe-3200" "TMe-2866"
#> [25] "TMe-2952" "TMe-3805" "TMe-2128" "TMe-1461" "TMe-3766" "TMe-1241"
#> [31] "TMe-86"  
#> 
#> [[3]]
#>  [1] "TMe-35"   "TMe-1897" "TMe-3638" "TMe-785"  "TMe-1593" "TMe-3147"
#>  [7] "TMe-3596" "TMe-946"  "TMe-3274" "TMe-3631" "TMe-141"  "TMe-1838"
#> [13] "TMe-234"  "TMe-2926" "TMe-1863" "TMe-2532" "TMe-3383" "TMe-200" 
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1434" "TMe-3527" "TMe-3255"
#>  [7] "TMe-3542" "TMe-650"  "TMe-2971" "TMe-2809" "TMe-663"  "TMe-1419"
#> [13] "TMe-1511" "TMe-3054" "TMe-3442" "TMe-1988" "TMe-241"  "TMe-2924"
#> [19] "TMe-2240" "TMe-2567" "TMe-3225" "TMe-3730" "TMe-3242" "TMe-550" 
#> [25] "TMe-1020" "TMe-225"  "TMe-1123" "TMe-3422" "TMe-1525" "TMe-3273"
#> [31] "TMe-1376" "TMe-2956" "TMe-3167" "TMe-388"  "TMe-3760"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2907" "TMe-3185" "TMe-277"  "TMe-1381"
#>  [7] "TMe-378"  "TMe-2953" "TMe-723"  "TMe-2589" "TMe-2531" "TMe-344" 
#> [13] "TMe-2853" "TMe-2121" "TMe-2041" "TMe-969"  "TMe-3034" "TMe-712" 
#> [19] "TMe-287"  "TMe-2790" "TMe-1500" "TMe-3356" "TMe-2355" "TMe-284" 
#> [25] "TMe-1004" "TMe-48"   "TMe-2820" "TMe-336"  "TMe-294"  "TMe-3628"
#> [31] "TMe-158"  "TMe-536"  "TMe-1458" "TMe-759"  "TMe-1924" "TMe-1963"
#> [37] "TMe-3311" "TMe-2296" "TMe-208"  "TMe-1357"
#> 
#> [[6]]
#>  [1] "TMe-725"  "TMe-1232" "TMe-1035" "TMe-3095" "TMe-1483" "TMe-693" 
#>  [7] "TMe-1428" "TMe-580"  "TMe-963"  "TMe-1875" "TMe-3549" "TMe-2510"
#> [13] "TMe-1353" "TMe-2791" "TMe-1124" "TMe-1816" "TMe-2064"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Pooled Brillouin diversity index")


# Mean Margalef's richness index
greedysel_first_mean_margalef <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_margalef,
                   metric = "mean", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_mean_margalef
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-3424" "TMe-778"  "TMe-2967"
#>  [7] "TMe-2965" "TMe-841"  "TMe-2010" "TMe-3030" "TMe-1140" "TMe-1915"
#> [13] "TMe-132"  "TMe-1226" "TMe-1360" "TMe-1564" "TMe-3262" "TMe-41"  
#> [19] "TMe-44"   "TMe-117"  "TMe-290"  "TMe-299"  "TMe-300"  "TMe-407" 
#> [25] "TMe-410"  "TMe-438"  "TMe-446"  "TMe-469"  "TMe-489"  "TMe-500" 
#> [31] "TMe-501" 
#> 
#> [[2]]
#>  [1] "TMe-377"  "TMe-1385" "TMe-477"  "TMe-2611" "TMe-3200" "TMe-509" 
#>  [7] "TMe-47"   "TMe-289"  "TMe-1907" "TMe-2033" "TMe-160"  "TMe-369" 
#> [13] "TMe-3466" "TMe-2304" "TMe-1241" "TMe-1461" "TMe-3222" "TMe-3237"
#> [19] "TMe-2"    "TMe-39"   "TMe-60"   "TMe-85"   "TMe-339"  "TMe-365" 
#> [25] "TMe-404"  "TMe-409"  "TMe-433"  "TMe-450"  "TMe-455"  "TMe-528" 
#> [31] "TMe-539" 
#> 
#> [[3]]
#>  [1] "TMe-434"  "TMe-773"  "TMe-3458" "TMe-3565" "TMe-405"  "TMe-946" 
#>  [7] "TMe-2897" "TMe-1897" "TMe-4"    "TMe-123"  "TMe-267"  "TMe-1868"
#> [13] "TMe-1993" "TMe-2592" "TMe-3046" "TMe-35"   "TMe-2926" "TMe-3596"
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1350" "TMe-552"  "TMe-1434"
#>  [7] "TMe-3527" "TMe-3054" "TMe-3255" "TMe-540"  "TMe-3428" "TMe-241" 
#> [13] "TMe-103"  "TMe-266"  "TMe-427"  "TMe-460"  "TMe-550"  "TMe-2809"
#> [19] "TMe-3273" "TMe-1020" "TMe-1511" "TMe-2567" "TMe-11"   "TMe-12"  
#> [25] "TMe-57"   "TMe-62"   "TMe-90"   "TMe-107"  "TMe-108"  "TMe-150" 
#> [31] "TMe-154"  "TMe-170"  "TMe-182"  "TMe-191"  "TMe-204" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2307" "TMe-3185" "TMe-197"  "TMe-723" 
#>  [7] "TMe-378"  "TMe-2953" "TMe-277"  "TMe-2853" "TMe-1004" "TMe-336" 
#> [13] "TMe-399"  "TMe-457"  "TMe-1904" "TMe-2041" "TMe-181"  "TMe-48"  
#> [19] "TMe-98"   "TMe-99"   "TMe-158"  "TMe-167"  "TMe-212"  "TMe-275" 
#> [25] "TMe-294"  "TMe-305"  "TMe-307"  "TMe-332"  "TMe-334"  "TMe-340" 
#> [31] "TMe-347"  "TMe-359"  "TMe-360"  "TMe-362"  "TMe-363"  "TMe-373" 
#> [37] "TMe-376"  "TMe-382"  "TMe-385"  "TMe-390" 
#> 
#> [[6]]
#>  [1] "TMe-856"  "TMe-3095" "TMe-693"  "TMe-1238" "TMe-1883" "TMe-2035"
#>  [7] "TMe-315"  "TMe-1858" "TMe-1875" "TMe-1124" "TMe-130"  "TMe-751" 
#> [13] "TMe-1775" "TMe-963"  "TMe-138"  "TMe-269"  "TMe-321" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Mean Margalef's richness index")


# Pooled Margalef's richness index
greedysel_first_sum_margalef <-
  select.diversity(data = data, names = "genotypes", group = "Cluster",
                   alloc = counts, qualitative = traits,
                   always.selected = mand_accns, div.fun = div_fun_margalef,
                   metric = "pooled", search = "greedy",
                   local.search = "first.improvement",max.iter = 3)
greedysel_first_sum_margalef
#> [[1]]
#>  [1] "TMe-3423" "TMe-2027" "TMe-867"  "TMe-677"  "TMe-3705" "TMe-1717"
#>  [7] "TMe-3030" "TMe-1425" "TMe-2010" "TMe-2069" "TMe-1140" "TMe-3396"
#> [13] "TMe-1226" "TMe-778"  "TMe-1360" "TMe-132"  "TMe-1922" "TMe-1915"
#> [19] "TMe-41"   "TMe-44"   "TMe-117"  "TMe-290"  "TMe-299"  "TMe-300" 
#> [25] "TMe-407"  "TMe-410"  "TMe-438"  "TMe-446"  "TMe-469"  "TMe-489" 
#> [31] "TMe-500" 
#> 
#> [[2]]
#>  [1] "TMe-3200" "TMe-2128" "TMe-2200" "TMe-3530" "TMe-289"  "TMe-3338"
#>  [7] "TMe-2033" "TMe-409"  "TMe-1907" "TMe-160"  "TMe-3466" "TMe-3237"
#> [13] "TMe-369"  "TMe-509"  "TMe-1241" "TMe-2304" "TMe-1461" "TMe-3222"
#> [19] "TMe-2"    "TMe-39"   "TMe-60"   "TMe-85"   "TMe-339"  "TMe-365" 
#> [25] "TMe-377"  "TMe-404"  "TMe-433"  "TMe-450"  "TMe-455"  "TMe-477" 
#> [31] "TMe-528" 
#> 
#> [[3]]
#>  [1] "TMe-3326" "TMe-1897" "TMe-785"  "TMe-3638" "TMe-773"  "TMe-946" 
#>  [7] "TMe-2977" "TMe-3100" "TMe-70"   "TMe-3046" "TMe-4"    "TMe-267" 
#> [13] "TMe-1868" "TMe-2926" "TMe-2592" "TMe-35"   "TMe-3596" "TMe-6"   
#> 
#> [[4]]
#>  [1] "TMe-34"   "TMe-801"  "TMe-698"  "TMe-1350" "TMe-552"  "TMe-1434"
#>  [7] "TMe-3255" "TMe-3527" "TMe-3054" "TMe-540"  "TMe-3428" "TMe-241" 
#> [13] "TMe-550"  "TMe-103"  "TMe-427"  "TMe-266"  "TMe-460"  "TMe-2809"
#> [19] "TMe-3273" "TMe-1020" "TMe-1511" "TMe-2567" "TMe-11"   "TMe-12"  
#> [25] "TMe-57"   "TMe-62"   "TMe-90"   "TMe-107"  "TMe-108"  "TMe-150" 
#> [31] "TMe-154"  "TMe-170"  "TMe-182"  "TMe-191"  "TMe-204" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-551"  "TMe-2307" "TMe-3185" "TMe-197"  "TMe-723" 
#>  [7] "TMe-378"  "TMe-2953" "TMe-277"  "TMe-2853" "TMe-1004" "TMe-336" 
#> [13] "TMe-399"  "TMe-457"  "TMe-1904" "TMe-2041" "TMe-181"  "TMe-48"  
#> [19] "TMe-98"   "TMe-99"   "TMe-158"  "TMe-167"  "TMe-212"  "TMe-275" 
#> [25] "TMe-294"  "TMe-305"  "TMe-307"  "TMe-332"  "TMe-334"  "TMe-340" 
#> [31] "TMe-347"  "TMe-359"  "TMe-360"  "TMe-362"  "TMe-363"  "TMe-373" 
#> [37] "TMe-376"  "TMe-382"  "TMe-385"  "TMe-390" 
#> 
#> [[6]]
#>  [1] "TMe-1883" "TMe-1518" "TMe-2963" "TMe-1509" "TMe-1858" "TMe-1174"
#>  [7] "TMe-696"  "TMe-1775" "TMe-3095" "TMe-315"  "TMe-751"  "TMe-1875"
#> [13] "TMe-236"  "TMe-130"  "TMe-856"  "TMe-138"  "TMe-269" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_richness,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Pooled Margalef's richness index")

```
