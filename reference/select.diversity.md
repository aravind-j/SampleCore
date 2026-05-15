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

set.seed(123)
cassava_EC_gp <- cassava_EC_gp[sample(1:nrow(cassava_EC_gp), 500), ]

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
  c("TMe-2018", "TMe-801", "TMe-3191", "TMe-1830", "TMe-1790")

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
#>  [1] "TMe-1830" "TMe-1922" "TMe-3112" "TMe-500"  "TMe-937"  "TMe-694" 
#>  [7] "TMe-3437" "TMe-2934" "TMe-815"  "TMe-2967" "TMe-3548" "TMe-2943"
#> [13] "TMe-1581" "TMe-3694" "TMe-1170" "TMe-2027" "TMe-1425" "TMe-3261"
#> [19] "TMe-3130" "TMe-44"   "TMe-2964" "TMe-3252" "TMe-1823" "TMe-300" 
#> [25] "TMe-1451" "TMe-1218" "TMe-3726" "TMe-2975" "TMe-865"  "TMe-1589"
#> [31] "TMe-2785"
#> 
#> [[2]]
#>  [1] "TMe-40"   "TMe-2951" "TMe-3530" "TMe-369"  "TMe-2611" "TMe-539" 
#>  [7] "TMe-2033" "TMe-1505" "TMe-890"  "TMe-3200" "TMe-3447" "TMe-674" 
#> [13] "TMe-74"   "TMe-196"  "TMe-3239" "TMe-2952" "TMe-2211" "TMe-3284"
#> [19] "TMe-3101" "TMe-960"  "TMe-2715" "TMe-2258" "TMe-455"  "TMe-1732"
#> [25] "TMe-409"  "TMe-2329" "TMe-2412" "TMe-3800" "TMe-1831" "TMe-3547"
#> [31] "TMe-2021"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3569" "TMe-3133" "TMe-261"  "TMe-3750" "TMe-1965"
#>  [7] "TMe-3397" "TMe-14"   "TMe-3556" "TMe-3336" "TMe-1804" "TMe-425" 
#> [13] "TMe-3620" "TMe-3638" "TMe-2394" "TMe-3644" "TMe-3715" "TMe-1725"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2788" "TMe-1016" "TMe-2958" "TMe-3558"
#>  [7] "TMe-27"   "TMe-25"   "TMe-460"  "TMe-3459" "TMe-1434" "TMe-3409"
#> [13] "TMe-812"  "TMe-1020" "TMe-1330" "TMe-3269" "TMe-2567" "TMe-2367"
#> [19] "TMe-191"  "TMe-2059" "TMe-2458" "TMe-2776" "TMe-1700" "TMe-619" 
#> [25] "TMe-2039" "TMe-3581" "TMe-622"  "TMe-241"  "TMe-1155" "TMe-1179"
#> [31] "TMe-3198" "TMe-3428" "TMe-3000" "TMe-368"  "TMe-2839"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-823"  "TMe-1037" "TMe-745"  "TMe-1559" "TMe-312" 
#>  [7] "TMe-307"  "TMe-348"  "TMe-574"  "TMe-682"  "TMe-1388" "TMe-2907"
#> [13] "TMe-1760" "TMe-712"  "TMe-362"  "TMe-3329" "TMe-419"  "TMe-2853"
#> [19] "TMe-730"  "TMe-769"  "TMe-1500" "TMe-688"  "TMe-920"  "TMe-1196"
#> [25] "TMe-2425" "TMe-600"  "TMe-1004" "TMe-1680" "TMe-629"  "TMe-603" 
#> [31] "TMe-2271" "TMe-487"  "TMe-1401" "TMe-588"  "TMe-423"  "TMe-1003"
#> [37] "TMe-360"  "TMe-1979" "TMe-816"  "TMe-1414"
#> 
#> [[6]]
#>  [1] "TMe-1445" "TMe-1413" "TMe-310"  "TMe-661"  "TMe-1945" "TMe-531" 
#>  [7] "TMe-1124" "TMe-751"  "TMe-2791" "TMe-1816" "TMe-2818" "TMe-693" 
#> [13] "TMe-1035" "TMe-1566" "TMe-1217" "TMe-620"  "TMe-598" 
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
#>  [1] "TMe-1830" "TMe-2237" "TMe-2944" "TMe-1091" "TMe-778"  "TMe-1425"
#>  [7] "TMe-3252" "TMe-3104" "TMe-3539" "TMe-3514" "TMe-3353" "TMe-882" 
#> [13] "TMe-44"   "TMe-910"  "TMe-1589" "TMe-2152" "TMe-1190" "TMe-937" 
#> [19] "TMe-3051" "TMe-3419" "TMe-694"  "TMe-2810" "TMe-1451" "TMe-3112"
#> [25] "TMe-3110" "TMe-566"  "TMe-3262" "TMe-2027" "TMe-2785" "TMe-3465"
#> [31] "TMe-3719"
#> 
#> [[2]]
#>  [1] "TMe-2568" "TMe-2033" "TMe-1474" "TMe-2329" "TMe-2021" "TMe-1323"
#>  [7] "TMe-890"  "TMe-960"  "TMe-74"   "TMe-539"  "TMe-3239" "TMe-1385"
#> [13] "TMe-3366" "TMe-2258" "TMe-2611" "TMe-3101" "TMe-2952" "TMe-369" 
#> [19] "TMe-3557" "TMe-40"   "TMe-171"  "TMe-409"  "TMe-1505" "TMe-796" 
#> [25] "TMe-3766" "TMe-196"  "TMe-3200" "TMe-251"  "TMe-2903" "TMe-3547"
#> [31] "TMe-2757"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-1965" "TMe-3556" "TMe-3569" "TMe-1738" "TMe-2086"
#>  [7] "TMe-161"  "TMe-3407" "TMe-234"  "TMe-3397" "TMe-381"  "TMe-3148"
#> [13] "TMe-14"   "TMe-155"  "TMe-1230" "TMe-123"  "TMe-1804" "TMe-3643"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2039" "TMe-191"  "TMe-1377" "TMe-18"  
#>  [7] "TMe-63"   "TMe-3189" "TMe-3269" "TMe-2958" "TMe-1020" "TMe-78"  
#> [13] "TMe-3327" "TMe-3072" "TMe-1144" "TMe-3248" "TMe-27"   "TMe-2839"
#> [19] "TMe-2788" "TMe-761"  "TMe-428"  "TMe-368"  "TMe-526"  "TMe-3390"
#> [25] "TMe-3378" "TMe-699"  "TMe-3459" "TMe-3527" "TMe-1434" "TMe-186" 
#> [31] "TMe-427"  "TMe-2841" "TMe-1297" "TMe-3581" "TMe-956" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-823"  "TMe-1004" "TMe-795"  "TMe-532"  "TMe-723" 
#>  [7] "TMe-1375" "TMe-1414" "TMe-1159" "TMe-712"  "TMe-600"  "TMe-1220"
#> [13] "TMe-1694" "TMe-803"  "TMe-627"  "TMe-1500" "TMe-755"  "TMe-892" 
#> [19] "TMe-543"  "TMe-707"  "TMe-816"  "TMe-1401" "TMe-2271" "TMe-877" 
#> [25] "TMe-1199" "TMe-2138" "TMe-2124" "TMe-645"  "TMe-167"  "TMe-2003"
#> [31] "TMe-1427" "TMe-769"  "TMe-1979" "TMe-920"  "TMe-623"  "TMe-423" 
#> [37] "TMe-588"  "TMe-98"   "TMe-332"  "TMe-312" 
#> 
#> [[6]]
#>  [1] "TMe-1744" "TMe-2818" "TMe-598"  "TMe-661"  "TMe-1124" "TMe-854" 
#>  [7] "TMe-1756" "TMe-1945" "TMe-1592" "TMe-2791" "TMe-1661" "TMe-1676"
#> [13] "TMe-1900" "TMe-693"  "TMe-1392" "TMe-1302" "TMe-968" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_sum_richness,
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
#>  [1] "TMe-1830" "TMe-3051" "TMe-778"  "TMe-1190" "TMe-1425" "TMe-3236"
#>  [7] "TMe-306"  "TMe-3262" "TMe-1869" "TMe-3481" "TMe-300"  "TMe-3685"
#> [13] "TMe-569"  "TMe-1218" "TMe-1717" "TMe-815"  "TMe-3726" "TMe-3719"
#> [19] "TMe-1930" "TMe-28"   "TMe-1451" "TMe-3353" "TMe-2934" "TMe-2066"
#> [25] "TMe-2964" "TMe-694"  "TMe-2027" "TMe-3514" "TMe-865"  "TMe-3111"
#> [31] "TMe-41"  
#> 
#> [[2]]
#>  [1] "TMe-409"  "TMe-369"  "TMe-2568" "TMe-2757" "TMe-1732" "TMe-3766"
#>  [7] "TMe-3547" "TMe-2000" "TMe-171"  "TMe-1754" "TMe-3101" "TMe-2257"
#> [13] "TMe-2412" "TMe-2903" "TMe-3200" "TMe-1385" "TMe-1668" "TMe-3530"
#> [19] "TMe-74"   "TMe-1474" "TMe-2952" "TMe-2997" "TMe-251"  "TMe-539" 
#> [25] "TMe-2033" "TMe-2352" "TMe-3093" "TMe-681"  "TMe-674"  "TMe-339" 
#> [31] "TMe-2021"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-2161" "TMe-3336" "TMe-3715" "TMe-14"   "TMe-2502"
#>  [7] "TMe-1725" "TMe-2086" "TMe-1819" "TMe-234"  "TMe-3631" "TMe-1176"
#> [13] "TMe-3207" "TMe-425"  "TMe-3133" "TMe-70"   "TMe-3620" "TMe-3638"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-280"  "TMe-3189" "TMe-3562" "TMe-3428"
#>  [7] "TMe-1456" "TMe-3558" "TMe-2947" "TMe-1434" "TMe-3242" "TMe-698" 
#> [13] "TMe-3585" "TMe-3591" "TMe-2928" "TMe-3000" "TMe-386"  "TMe-3527"
#> [19] "TMe-812"  "TMe-279"  "TMe-875"  "TMe-3273" "TMe-622"  "TMe-460" 
#> [25] "TMe-897"  "TMe-761"  "TMe-2839" "TMe-3511" "TMe-427"  "TMe-372" 
#> [31] "TMe-241"  "TMe-1150" "TMe-1246" "TMe-3390" "TMe-1377"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-574"  "TMe-2753" "TMe-1534" "TMe-2192" "TMe-2138"
#>  [7] "TMe-1880" "TMe-1401" "TMe-2853" "TMe-1500" "TMe-870"  "TMe-1037"
#> [13] "TMe-623"  "TMe-1131" "TMe-1388" "TMe-247"  "TMe-2290" "TMe-712" 
#> [19] "TMe-1159" "TMe-2907" "TMe-2003" "TMe-1427" "TMe-682"  "TMe-730" 
#> [25] "TMe-1541" "TMe-1453" "TMe-406"  "TMe-1220" "TMe-667"  "TMe-307" 
#> [31] "TMe-1414" "TMe-167"  "TMe-1160" "TMe-344"  "TMe-600"  "TMe-2355"
#> [37] "TMe-723"  "TMe-892"  "TMe-1234" "TMe-769" 
#> 
#> [[6]]
#>  [1] "TMe-598"  "TMe-1676" "TMe-1403" "TMe-531"  "TMe-1661" "TMe-854" 
#>  [7] "TMe-2983" "TMe-1945" "TMe-310"  "TMe-2196" "TMe-1900" "TMe-631" 
#> [13] "TMe-693"  "TMe-1124" "TMe-1816" "TMe-2791" "TMe-620" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_shannon,
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
#>  [1] "TMe-1830" "TMe-3337" "TMe-3262" "TMe-2943" "TMe-1930" "TMe-3115"
#>  [7] "TMe-1190" "TMe-2964" "TMe-3437" "TMe-3236" "TMe-2785" "TMe-1086"
#> [13] "TMe-1425" "TMe-937"  "TMe-2934" "TMe-1451" "TMe-882"  "TMe-44"  
#> [19] "TMe-1960" "TMe-952"  "TMe-41"   "TMe-3514" "TMe-2810" "TMe-3623"
#> [25] "TMe-3726" "TMe-2103" "TMe-306"  "TMe-815"  "TMe-28"   "TMe-2967"
#> [31] "TMe-1117"
#> 
#> [[2]]
#>  [1] "TMe-960"  "TMe-2951" "TMe-2952" "TMe-3547" "TMe-3800" "TMe-3557"
#>  [7] "TMe-3101" "TMe-2352" "TMe-2211" "TMe-1668" "TMe-251"  "TMe-2000"
#> [13] "TMe-1323" "TMe-74"   "TMe-3530" "TMe-2329" "TMe-3239" "TMe-1754"
#> [19] "TMe-2033" "TMe-3200" "TMe-2257" "TMe-1385" "TMe-369"  "TMe-796" 
#> [25] "TMe-2021" "TMe-3447" "TMe-2568" "TMe-1831" "TMe-2715" "TMe-2412"
#> [31] "TMe-890" 
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-1230" "TMe-70"   "TMe-2394" "TMe-1863" "TMe-2502"
#>  [7] "TMe-1939" "TMe-3397" "TMe-3644" "TMe-234"  "TMe-3643" "TMe-830" 
#> [13] "TMe-3148" "TMe-3230" "TMe-1804" "TMe-3638" "TMe-425"  "TMe-3556"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-3108" "TMe-1139" "TMe-1988" "TMe-12"  
#>  [7] "TMe-1123" "TMe-1150" "TMe-266"  "TMe-619"  "TMe-1579" "TMe-3428"
#> [13] "TMe-699"  "TMe-2841" "TMe-170"  "TMe-2924" "TMe-1148" "TMe-3198"
#> [19] "TMe-3542" "TMe-1776" "TMe-3265" "TMe-318"  "TMe-2971" "TMe-1297"
#> [25] "TMe-1987" "TMe-259"  "TMe-3527" "TMe-1350" "TMe-2755" "TMe-897" 
#> [31] "TMe-594"  "TMe-3242" "TMe-3558" "TMe-2240" "TMe-2928"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-1037" "TMe-1375" "TMe-769"  "TMe-1440" "TMe-863" 
#>  [7] "TMe-439"  "TMe-1003" "TMe-2003" "TMe-803"  "TMe-167"  "TMe-98"  
#> [13] "TMe-1295" "TMe-2290" "TMe-920"  "TMe-712"  "TMe-1196" "TMe-2425"
#> [19] "TMe-1159" "TMe-688"  "TMe-627"  "TMe-2138" "TMe-823"  "TMe-1979"
#> [25] "TMe-1427" "TMe-2750" "TMe-1293" "TMe-1934" "TMe-795"  "TMe-419" 
#> [31] "TMe-2907" "TMe-1414" "TMe-2855" "TMe-682"  "TMe-600"  "TMe-1299"
#> [37] "TMe-585"  "TMe-1788" "TMe-487"  "TMe-954" 
#> 
#> [[6]]
#>  [1] "TMe-1832" "TMe-1035" "TMe-310"  "TMe-2791" "TMe-3177" "TMe-505" 
#>  [7] "TMe-1124" "TMe-1566" "TMe-1445" "TMe-936"  "TMe-1608" "TMe-751" 
#> [13] "TMe-1302" "TMe-693"  "TMe-1661" "TMe-1503" "TMe-2543"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_sum_shannon,
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
#>  [1] "TMe-1830" "TMe-28"   "TMe-1425" "TMe-2943" "TMe-815"  "TMe-3719"
#>  [7] "TMe-1960" "TMe-3104" "TMe-3726" "TMe-569"  "TMe-2103" "TMe-3351"
#> [13] "TMe-3501" "TMe-1086" "TMe-3252" "TMe-300"  "TMe-1451" "TMe-2237"
#> [19] "TMe-1532" "TMe-2513" "TMe-3437" "TMe-3051" "TMe-500"  "TMe-1190"
#> [25] "TMe-3353" "TMe-306"  "TMe-2967" "TMe-3685" "TMe-41"   "TMe-1589"
#> [31] "TMe-1581"
#> 
#> [[2]]
#>  [1] "TMe-1323" "TMe-1698" "TMe-3239" "TMe-960"  "TMe-1668" "TMe-369" 
#>  [7] "TMe-171"  "TMe-3805" "TMe-1831" "TMe-2329" "TMe-3766" "TMe-1754"
#> [13] "TMe-1732" "TMe-2995" "TMe-3200" "TMe-2033" "TMe-3530" "TMe-2952"
#> [19] "TMe-2903" "TMe-2715" "TMe-3495" "TMe-2352" "TMe-40"   "TMe-890" 
#> [25] "TMe-2000" "TMe-2611" "TMe-3366" "TMe-196"  "TMe-3284" "TMe-3101"
#> [31] "TMe-3258"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-1198" "TMe-381"  "TMe-3638" "TMe-3230" "TMe-1863"
#>  [7] "TMe-3148" "TMe-1230" "TMe-3407" "TMe-2161" "TMe-425"  "TMe-3643"
#> [13] "TMe-234"  "TMe-3750" "TMe-70"   "TMe-2502" "TMe-1787" "TMe-3556"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-3390" "TMe-2776" "TMe-1297" "TMe-1988"
#>  [7] "TMe-3000" "TMe-3409" "TMe-3542" "TMe-2247" "TMe-3558" "TMe-259" 
#> [13] "TMe-761"  "TMe-3591" "TMe-1150" "TMe-241"  "TMe-1148" "TMe-698" 
#> [19] "TMe-3108" "TMe-2755" "TMe-1179" "TMe-186"  "TMe-1336" "TMe-516" 
#> [25] "TMe-1456" "TMe-1350" "TMe-3198" "TMe-3585" "TMe-65"   "TMe-279" 
#> [31] "TMe-3538" "TMe-1027" "TMe-1434" "TMe-875"  "TMe-3273"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-730"  "TMe-2425" "TMe-2192" "TMe-288"  "TMe-362" 
#>  [7] "TMe-723"  "TMe-344"  "TMe-1290" "TMe-1388" "TMe-1694" "TMe-332" 
#> [13] "TMe-2151" "TMe-712"  "TMe-629"  "TMe-645"  "TMe-823"  "TMe-2355"
#> [19] "TMe-1234" "TMe-256"  "TMe-419"  "TMe-1188" "TMe-2590" "TMe-1003"
#> [25] "TMe-1500" "TMe-2853" "TMe-307"  "TMe-954"  "TMe-627"  "TMe-1622"
#> [31] "TMe-667"  "TMe-863"  "TMe-1295" "TMe-323"  "TMe-2855" "TMe-224" 
#> [37] "TMe-98"   "TMe-803"  "TMe-1196" "TMe-2271"
#> 
#> [[6]]
#>  [1] "TMe-1566" "TMe-854"  "TMe-310"  "TMe-2818" "TMe-1661" "TMe-2791"
#>  [7] "TMe-2983" "TMe-1608" "TMe-661"  "TMe-1302" "TMe-1445" "TMe-1503"
#> [13] "TMe-531"  "TMe-1124" "TMe-222"  "TMe-693"  "TMe-631" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_simpson,
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
#>  [1] "TMe-1830" "TMe-3514" "TMe-1218" "TMe-2967" "TMe-937"  "TMe-952" 
#>  [7] "TMe-3501" "TMe-1170" "TMe-2975" "TMe-3236" "TMe-3252" "TMe-1425"
#> [13] "TMe-2152" "TMe-3623" "TMe-1190" "TMe-44"   "TMe-2103" "TMe-1096"
#> [19] "TMe-569"  "TMe-566"  "TMe-2785" "TMe-815"  "TMe-3694" "TMe-41"  
#> [25] "TMe-2027" "TMe-3548" "TMe-3112" "TMe-3142" "TMe-3719" "TMe-1973"
#> [31] "TMe-3419"
#> 
#> [[2]]
#>  [1] "TMe-3366" "TMe-3530" "TMe-2352" "TMe-3101" "TMe-369"  "TMe-1732"
#>  [7] "TMe-3200" "TMe-1505" "TMe-2568" "TMe-1754" "TMe-890"  "TMe-3093"
#> [13] "TMe-2329" "TMe-2211" "TMe-3239" "TMe-2715" "TMe-40"   "TMe-1385"
#> [19] "TMe-3805" "TMe-1474" "TMe-2611" "TMe-2412" "TMe-1668" "TMe-2258"
#> [25] "TMe-3447" "TMe-674"  "TMe-681"  "TMe-2951" "TMe-539"  "TMe-2033"
#> [31] "TMe-3800"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3100" "TMe-70"   "TMe-3750" "TMe-14"   "TMe-1863"
#>  [7] "TMe-2748" "TMe-1230" "TMe-261"  "TMe-3620" "TMe-1819" "TMe-3397"
#> [13] "TMe-3088" "TMe-2086" "TMe-381"  "TMe-3556" "TMe-2968" "TMe-161" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-3269" "TMe-25"   "TMe-1376" "TMe-761" 
#>  [7] "TMe-2956" "TMe-2458" "TMe-63"   "TMe-3257" "TMe-3189" "TMe-3327"
#> [13] "TMe-78"   "TMe-280"  "TMe-353"  "TMe-1148" "TMe-2841" "TMe-387" 
#> [19] "TMe-386"  "TMe-1123" "TMe-65"   "TMe-368"  "TMe-1179" "TMe-698" 
#> [25] "TMe-812"  "TMe-3542" "TMe-1700" "TMe-2958" "TMe-2971" "TMe-57"  
#> [31] "TMe-3242" "TMe-259"  "TMe-956"  "TMe-3511" "TMe-1020"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-1012" "TMe-863"  "TMe-487"  "TMe-224"  "TMe-323" 
#>  [7] "TMe-1880" "TMe-1299" "TMe-1265" "TMe-2355" "TMe-627"  "TMe-2192"
#> [13] "TMe-3329" "TMe-816"  "TMe-333"  "TMe-603"  "TMe-2151" "TMe-2271"
#> [19] "TMe-1293" "TMe-532"  "TMe-344"  "TMe-1453" "TMe-1131" "TMe-574" 
#> [25] "TMe-439"  "TMe-1788" "TMe-920"  "TMe-2016" "TMe-1680" "TMe-688" 
#> [31] "TMe-755"  "TMe-423"  "TMe-2907" "TMe-1196" "TMe-256"  "TMe-1694"
#> [37] "TMe-730"  "TMe-892"  "TMe-2590" "TMe-2855"
#> 
#> [[6]]
#>  [1] "TMe-310"  "TMe-1816" "TMe-620"  "TMe-1124" "TMe-1035" "TMe-2543"
#>  [7] "TMe-2196" "TMe-1403" "TMe-598"  "TMe-1413" "TMe-1608" "TMe-1592"
#> [13] "TMe-1566" "TMe-1646" "TMe-1302" "TMe-1661" "TMe-2791"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_sum_simpson,
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
#>  [1] "TMe-1830" "TMe-1922" "TMe-2237" "TMe-3419" "TMe-3539" "TMe-1869"
#>  [7] "TMe-469"  "TMe-3252" "TMe-3726" "TMe-41"   "TMe-3142" "TMe-1960"
#> [13] "TMe-3111" "TMe-1091" "TMe-882"  "TMe-778"  "TMe-1117" "TMe-3262"
#> [19] "TMe-3236" "TMe-28"   "TMe-1717" "TMe-3465" "TMe-1451" "TMe-44"  
#> [25] "TMe-2967" "TMe-306"  "TMe-1190" "TMe-2943" "TMe-2066" "TMe-3115"
#> [31] "TMe-3337"
#> 
#> [[2]]
#>  [1] "TMe-1698" "TMe-539"  "TMe-409"  "TMe-339"  "TMe-2757" "TMe-3766"
#>  [7] "TMe-1831" "TMe-2568" "TMe-2352" "TMe-3258" "TMe-2329" "TMe-3200"
#> [13] "TMe-2611" "TMe-3101" "TMe-3530" "TMe-251"  "TMe-2715" "TMe-2951"
#> [19] "TMe-369"  "TMe-3495" "TMe-2903" "TMe-1385" "TMe-2000" "TMe-2258"
#> [25] "TMe-2033" "TMe-1754" "TMe-3805" "TMe-2952" "TMe-196"  "TMe-1732"
#> [31] "TMe-674" 
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-425"  "TMe-3336" "TMe-3750" "TMe-3556" "TMe-155" 
#>  [7] "TMe-2374" "TMe-2086" "TMe-3133" "TMe-1863" "TMe-234"  "TMe-381" 
#> [13] "TMe-2169" "TMe-261"  "TMe-830"  "TMe-3230" "TMe-1230" "TMe-14"  
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-594"  "TMe-3602" "TMe-1144" "TMe-1434"
#>  [7] "TMe-1139" "TMe-3103" "TMe-3327" "TMe-3255" "TMe-3459" "TMe-3189"
#> [13] "TMe-1579" "TMe-2247" "TMe-812"  "TMe-266"  "TMe-2240" "TMe-1988"
#> [19] "TMe-3257" "TMe-372"  "TMe-1456" "TMe-2567" "TMe-1155" "TMe-3390"
#> [25] "TMe-3144" "TMe-2788" "TMe-3273" "TMe-1246" "TMe-2841" "TMe-3248"
#> [31] "TMe-25"   "TMe-1700" "TMe-2059" "TMe-3537" "TMe-2458"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-1012" "TMe-247"  "TMe-245"  "TMe-2003" "TMe-2355"
#>  [7] "TMe-2853" "TMe-224"  "TMe-406"  "TMe-730"  "TMe-892"  "TMe-877" 
#> [13] "TMe-360"  "TMe-2441" "TMe-574"  "TMe-1500" "TMe-920"  "TMe-667" 
#> [19] "TMe-755"  "TMe-870"  "TMe-1196" "TMe-1534" "TMe-2907" "TMe-2590"
#> [25] "TMe-1788" "TMe-1934" "TMe-585"  "TMe-362"  "TMe-1160" "TMe-623" 
#> [31] "TMe-2290" "TMe-2425" "TMe-1220" "TMe-1159" "TMe-600"  "TMe-2057"
#> [37] "TMe-1760" "TMe-823"  "TMe-307"  "TMe-167" 
#> 
#> [[6]]
#>  [1] "TMe-1392" "TMe-505"  "TMe-310"  "TMe-1302" "TMe-1900" "TMe-2791"
#>  [7] "TMe-2543" "TMe-1992" "TMe-2818" "TMe-1816" "TMe-2983" "TMe-1076"
#> [13] "TMe-854"  "TMe-3177" "TMe-1566" "TMe-1832" "TMe-1035"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_mcintosh,
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
#>  [1] "TMe-1830" "TMe-1086" "TMe-469"  "TMe-3694" "TMe-3465" "TMe-3142"
#>  [7] "TMe-2785" "TMe-2934" "TMe-41"   "TMe-2152" "TMe-3623" "TMe-1451"
#> [13] "TMe-3685" "TMe-500"  "TMe-3236" "TMe-3514" "TMe-1589" "TMe-2943"
#> [19] "TMe-815"  "TMe-865"  "TMe-1717" "TMe-3539" "TMe-1869" "TMe-3111"
#> [25] "TMe-1170" "TMe-778"  "TMe-3437" "TMe-3419" "TMe-569"  "TMe-1581"
#> [31] "TMe-3726"
#> 
#> [[2]]
#>  [1] "TMe-2258" "TMe-681"  "TMe-171"  "TMe-369"  "TMe-890"  "TMe-1323"
#>  [7] "TMe-3766" "TMe-1385" "TMe-74"   "TMe-2952" "TMe-3547" "TMe-1732"
#> [13] "TMe-3805" "TMe-455"  "TMe-796"  "TMe-539"  "TMe-251"  "TMe-674" 
#> [19] "TMe-3366" "TMe-2257" "TMe-2000" "TMe-3258" "TMe-2211" "TMe-3284"
#> [25] "TMe-40"   "TMe-1668" "TMe-3200" "TMe-2568" "TMe-3101" "TMe-2329"
#> [31] "TMe-2997"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-123"  "TMe-3335" "TMe-3631" "TMe-1804" "TMe-234" 
#>  [7] "TMe-3336" "TMe-2968" "TMe-3397" "TMe-3148" "TMe-3620" "TMe-1863"
#> [13] "TMe-3715" "TMe-14"   "TMe-155"  "TMe-2733" "TMe-3407" "TMe-3750"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-1776" "TMe-3378" "TMe-2367" "TMe-3255"
#>  [7] "TMe-3581" "TMe-516"  "TMe-3242" "TMe-3198" "TMe-1336" "TMe-2059"
#> [13] "TMe-1297" "TMe-1996" "TMe-3428" "TMe-3218" "TMe-698"  "TMe-1016"
#> [19] "TMe-1377" "TMe-1434" "TMe-1350" "TMe-2240" "TMe-3000" "TMe-525" 
#> [25] "TMe-3537" "TMe-27"   "TMe-2841" "TMe-1988" "TMe-3189" "TMe-956" 
#> [31] "TMe-3602" "TMe-1456" "TMe-3591" "TMe-3327" "TMe-2567"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-362"  "TMe-730"  "TMe-487"  "TMe-98"   "TMe-2590"
#>  [7] "TMe-863"  "TMe-256"  "TMe-2057" "TMe-3329" "TMe-823"  "TMe-997" 
#> [13] "TMe-1559" "TMe-1196" "TMe-755"  "TMe-1534" "TMe-1290" "TMe-629" 
#> [19] "TMe-1199" "TMe-954"  "TMe-419"  "TMe-1880" "TMe-712"  "TMe-2907"
#> [25] "TMe-1680" "TMe-1401" "TMe-2003" "TMe-333"  "TMe-1694" "TMe-167" 
#> [31] "TMe-647"  "TMe-2213" "TMe-2425" "TMe-600"  "TMe-307"  "TMe-1188"
#> [37] "TMe-1234" "TMe-816"  "TMe-1541" "TMe-423" 
#> 
#> [[6]]
#>  [1] "TMe-1900" "TMe-1566" "TMe-1302" "TMe-2818" "TMe-505"  "TMe-1756"
#>  [7] "TMe-1445" "TMe-1661" "TMe-2983" "TMe-661"  "TMe-2543" "TMe-598" 
#> [13] "TMe-693"  "TMe-1403" "TMe-1062" "TMe-2791" "TMe-1592"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_sum_mcintosh,
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
#>  [1] "TMe-1830" "TMe-3437" "TMe-1091" "TMe-1922" "TMe-300"  "TMe-2934"
#>  [7] "TMe-28"   "TMe-3481" "TMe-1425" "TMe-1096" "TMe-2944" "TMe-2785"
#> [13] "TMe-3252" "TMe-778"  "TMe-569"  "TMe-3236" "TMe-3353" "TMe-1581"
#> [19] "TMe-1930" "TMe-1960" "TMe-3261" "TMe-3104" "TMe-500"  "TMe-1190"
#> [25] "TMe-1914" "TMe-44"   "TMe-41"   "TMe-1117" "TMe-2066" "TMe-3501"
#> [31] "TMe-694" 
#> 
#> [[2]]
#>  [1] "TMe-74"   "TMe-2211" "TMe-3258" "TMe-2033" "TMe-674"  "TMe-2568"
#>  [7] "TMe-2257" "TMe-2903" "TMe-40"   "TMe-3366" "TMe-1474" "TMe-1732"
#> [13] "TMe-2000" "TMe-1385" "TMe-796"  "TMe-3800" "TMe-3200" "TMe-890" 
#> [19] "TMe-2329" "TMe-409"  "TMe-171"  "TMe-1754" "TMe-3447" "TMe-681" 
#> [25] "TMe-2021" "TMe-3766" "TMe-2258" "TMe-369"  "TMe-455"  "TMe-196" 
#> [31] "TMe-3101"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-1198" "TMe-3556" "TMe-234"  "TMe-3631" "TMe-1819"
#>  [7] "TMe-2086" "TMe-1230" "TMe-70"   "TMe-3007" "TMe-1863" "TMe-3643"
#> [13] "TMe-2394" "TMe-381"  "TMe-1939" "TMe-3335" "TMe-425"  "TMe-3750"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-3068" "TMe-1246" "TMe-3000" "TMe-1776"
#>  [7] "TMe-3581" "TMe-27"   "TMe-737"  "TMe-1579" "TMe-3072" "TMe-761" 
#> [13] "TMe-1376" "TMe-2567" "TMe-1988" "TMe-2946" "TMe-3428" "TMe-2924"
#> [19] "TMe-25"   "TMe-3250" "TMe-2039" "TMe-3103" "TMe-3257" "TMe-2059"
#> [25] "TMe-3232" "TMe-266"  "TMe-619"  "TMe-2318" "TMe-1144" "TMe-353" 
#> [31] "TMe-1027" "TMe-1020" "TMe-12"   "TMe-3409" "TMe-2788"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-707"  "TMe-1188" "TMe-1160" "TMe-2003" "TMe-439" 
#>  [7] "TMe-1934" "TMe-2271" "TMe-920"  "TMe-167"  "TMe-600"  "TMe-423" 
#> [13] "TMe-730"  "TMe-1293" "TMe-877"  "TMe-1375" "TMe-307"  "TMe-1220"
#> [19] "TMe-803"  "TMe-795"  "TMe-863"  "TMe-1265" "TMe-2355" "TMe-2192"
#> [25] "TMe-1453" "TMe-651"  "TMe-823"  "TMe-627"  "TMe-3329" "TMe-1440"
#> [31] "TMe-997"  "TMe-1159" "TMe-348"  "TMe-532"  "TMe-1234" "TMe-1401"
#> [37] "TMe-1199" "TMe-723"  "TMe-419"  "TMe-543" 
#> 
#> [[6]]
#>  [1] "TMe-661"  "TMe-1566" "TMe-693"  "TMe-222"  "TMe-1413" "TMe-1608"
#>  [7] "TMe-2818" "TMe-2983" "TMe-1445" "TMe-1403" "TMe-1503" "TMe-1392"
#> [13] "TMe-2791" "TMe-1816" "TMe-531"  "TMe-1992" "TMe-1945"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_brillouin,
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
#>  [1] "TMe-1830" "TMe-1091" "TMe-865"  "TMe-2943" "TMe-3351" "TMe-3115"
#>  [7] "TMe-1589" "TMe-3262" "TMe-1532" "TMe-1425" "TMe-888"  "TMe-1096"
#> [13] "TMe-2934" "TMe-3623" "TMe-3465" "TMe-44"   "TMe-3437" "TMe-3337"
#> [19] "TMe-3236" "TMe-1451" "TMe-28"   "TMe-952"  "TMe-1823" "TMe-2967"
#> [25] "TMe-1117" "TMe-3252" "TMe-2453" "TMe-2944" "TMe-3481" "TMe-500" 
#> [31] "TMe-41"  
#> 
#> [[2]]
#>  [1] "TMe-2952" "TMe-3093" "TMe-2412" "TMe-2000" "TMe-171"  "TMe-2903"
#>  [7] "TMe-2568" "TMe-681"  "TMe-369"  "TMe-1698" "TMe-2715" "TMe-3805"
#> [13] "TMe-1732" "TMe-960"  "TMe-674"  "TMe-1505" "TMe-3101" "TMe-3200"
#> [19] "TMe-3557" "TMe-539"  "TMe-1385" "TMe-1831" "TMe-2033" "TMe-2352"
#> [25] "TMe-3766" "TMe-2021" "TMe-3366" "TMe-2257" "TMe-3547" "TMe-2995"
#> [31] "TMe-1323"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-1804" "TMe-3407" "TMe-3230" "TMe-3336" "TMe-2169"
#>  [7] "TMe-3631" "TMe-1230" "TMe-64"   "TMe-1863" "TMe-2748" "TMe-2086"
#> [13] "TMe-2270" "TMe-381"  "TMe-261"  "TMe-1198" "TMe-155"  "TMe-1965"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-526"  "TMe-1776" "TMe-1336" "TMe-2247"
#>  [7] "TMe-3327" "TMe-3409" "TMe-956"  "TMe-1123" "TMe-3144" "TMe-266" 
#> [13] "TMe-3248" "TMe-2956" "TMe-1297" "TMe-1434" "TMe-3537" "TMe-3542"
#> [19] "TMe-1377" "TMe-3538" "TMe-57"   "TMe-619"  "TMe-428"  "TMe-1988"
#> [25] "TMe-2788" "TMe-897"  "TMe-2841" "TMe-78"   "TMe-241"  "TMe-3511"
#> [31] "TMe-3265" "TMe-2776" "TMe-3108" "TMe-698"  "TMe-3527"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-344"  "TMe-629"  "TMe-712"  "TMe-600"  "TMe-1004"
#>  [7] "TMe-532"  "TMe-2016" "TMe-2355" "TMe-332"  "TMe-1760" "TMe-1159"
#> [13] "TMe-823"  "TMe-1299" "TMe-423"  "TMe-816"  "TMe-1295" "TMe-1440"
#> [19] "TMe-612"  "TMe-603"  "TMe-1559" "TMe-2192" "TMe-1934" "TMe-2124"
#> [25] "TMe-795"  "TMe-803"  "TMe-2003" "TMe-2271" "TMe-723"  "TMe-1500"
#> [31] "TMe-1234" "TMe-651"  "TMe-1694" "TMe-2853" "TMe-574"  "TMe-920" 
#> [37] "TMe-954"  "TMe-419"  "TMe-256"  "TMe-1680"
#> 
#> [[6]]
#>  [1] "TMe-310"  "TMe-1445" "TMe-620"  "TMe-693"  "TMe-1124" "TMe-1035"
#>  [7] "TMe-1392" "TMe-1403" "TMe-1816" "TMe-598"  "TMe-1239" "TMe-2983"
#> [13] "TMe-1062" "TMe-2543" "TMe-2818" "TMe-1832" "TMe-1992"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_sum_brillouin,
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
#>  [1] "TMe-1830" "TMe-3111" "TMe-2785" "TMe-910"  "TMe-3623" "TMe-569" 
#>  [7] "TMe-3419" "TMe-3398" "TMe-2943" "TMe-2066" "TMe-937"  "TMe-500" 
#> [13] "TMe-3051" "TMe-778"  "TMe-1922" "TMe-3252" "TMe-3337" "TMe-300" 
#> [19] "TMe-2103" "TMe-2513" "TMe-306"  "TMe-41"   "TMe-3539" "TMe-2152"
#> [25] "TMe-1218" "TMe-1096" "TMe-2810" "TMe-2237" "TMe-3112" "TMe-2934"
#> [31] "TMe-1425"
#> 
#> [[2]]
#>  [1] "TMe-1754" "TMe-2952" "TMe-3101" "TMe-3239" "TMe-539"  "TMe-2951"
#>  [7] "TMe-339"  "TMe-2211" "TMe-455"  "TMe-681"  "TMe-2033" "TMe-3093"
#> [13] "TMe-1505" "TMe-409"  "TMe-3258" "TMe-369"  "TMe-2000" "TMe-2329"
#> [19] "TMe-2352" "TMe-796"  "TMe-2257" "TMe-674"  "TMe-2903" "TMe-2611"
#> [25] "TMe-3200" "TMe-3800" "TMe-74"   "TMe-1385" "TMe-2995" "TMe-40"  
#> [31] "TMe-3547"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-830"  "TMe-381"  "TMe-234"  "TMe-14"   "TMe-1939"
#>  [7] "TMe-3715" "TMe-3335" "TMe-1863" "TMe-2394" "TMe-1176" "TMe-3100"
#> [13] "TMe-261"  "TMe-1804" "TMe-161"  "TMe-3638" "TMe-1819" "TMe-3569"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-318"  "TMe-2240" "TMe-3327" "TMe-3390"
#>  [7] "TMe-3255" "TMe-1148" "TMe-1123" "TMe-1350" "TMe-386"  "TMe-840" 
#> [13] "TMe-594"  "TMe-1027" "TMe-3602" "TMe-241"  "TMe-1020" "TMe-372" 
#> [19] "TMe-3409" "TMe-387"  "TMe-27"   "TMe-2788" "TMe-1987" "TMe-2458"
#> [25] "TMe-812"  "TMe-3269" "TMe-3196" "TMe-2375" "TMe-65"   "TMe-875" 
#> [31] "TMe-2039" "TMe-1434" "TMe-3537" "TMe-3758" "TMe-3108"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-2853" "TMe-247"  "TMe-1375" "TMe-651"  "TMe-1427"
#>  [7] "TMe-1788" "TMe-348"  "TMe-2441" "TMe-344"  "TMe-1220" "TMe-816" 
#> [13] "TMe-600"  "TMe-877"  "TMe-574"  "TMe-1004" "TMe-419"  "TMe-585" 
#> [19] "TMe-1265" "TMe-1541" "TMe-2750" "TMe-645"  "TMe-1500" "TMe-487" 
#> [25] "TMe-1037" "TMe-323"  "TMe-723"  "TMe-1003" "TMe-627"  "TMe-98"  
#> [31] "TMe-1934" "TMe-1199" "TMe-2907" "TMe-1760" "TMe-1453" "TMe-647" 
#> [37] "TMe-1622" "TMe-2590" "TMe-1979" "TMe-2192"
#> 
#> [[6]]
#>  [1] "TMe-1217" "TMe-531"  "TMe-1832" "TMe-1646" "TMe-693"  "TMe-1403"
#>  [7] "TMe-751"  "TMe-1945" "TMe-1302" "TMe-1661" "TMe-2983" "TMe-2791"
#> [13] "TMe-1503" "TMe-1239" "TMe-1035" "TMe-1445" "TMe-1124"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_mean_margalef,
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
#>  [1] "TMe-1830" "TMe-2944" "TMe-500"  "TMe-2453" "TMe-1451" "TMe-1218"
#>  [7] "TMe-469"  "TMe-1425" "TMe-3252" "TMe-3142" "TMe-3112" "TMe-28"  
#> [13] "TMe-3115" "TMe-2934" "TMe-1869" "TMe-2785" "TMe-3719" "TMe-1930"
#> [19] "TMe-44"   "TMe-3685" "TMe-2027" "TMe-3110" "TMe-778"  "TMe-3236"
#> [25] "TMe-888"  "TMe-1086" "TMe-3398" "TMe-3539" "TMe-2103" "TMe-1589"
#> [31] "TMe-3694"
#> 
#> [[2]]
#>  [1] "TMe-40"   "TMe-3101" "TMe-3200" "TMe-674"  "TMe-3239" "TMe-1831"
#>  [7] "TMe-2903" "TMe-74"   "TMe-2757" "TMe-3557" "TMe-3547" "TMe-960" 
#> [13] "TMe-2211" "TMe-339"  "TMe-2568" "TMe-2715" "TMe-369"  "TMe-455" 
#> [19] "TMe-3284" "TMe-251"  "TMe-1754" "TMe-2951" "TMe-3258" "TMe-409" 
#> [25] "TMe-3766" "TMe-2033" "TMe-539"  "TMe-2611" "TMe-3805" "TMe-1698"
#> [31] "TMe-3495"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3638" "TMe-3631" "TMe-1198" "TMe-2502" "TMe-234" 
#>  [7] "TMe-3230" "TMe-2394" "TMe-155"  "TMe-3715" "TMe-1863" "TMe-1819"
#> [13] "TMe-148"  "TMe-830"  "TMe-1230" "TMe-1738" "TMe-3556" "TMe-1804"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-3542" "TMe-3459" "TMe-699"  "TMe-1434"
#>  [7] "TMe-2755" "TMe-3428" "TMe-1350" "TMe-368"  "TMe-12"   "TMe-279" 
#> [13] "TMe-3068" "TMe-2971" "TMe-353"  "TMe-812"  "TMe-427"  "TMe-2839"
#> [19] "TMe-1167" "TMe-3585" "TMe-25"   "TMe-2788" "TMe-3378" "TMe-1179"
#> [25] "TMe-65"   "TMe-1336" "TMe-1150" "TMe-619"  "TMe-1376" "TMe-3144"
#> [31] "TMe-1700" "TMe-3232" "TMe-1988" "TMe-241"  "TMe-3527"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-323"  "TMe-2590" "TMe-98"   "TMe-1979" "TMe-627" 
#>  [7] "TMe-2138" "TMe-1131" "TMe-2057" "TMe-439"  "TMe-2192" "TMe-344" 
#> [13] "TMe-1880" "TMe-1196" "TMe-723"  "TMe-2853" "TMe-1622" "TMe-2855"
#> [19] "TMe-1159" "TMe-1414" "TMe-954"  "TMe-2907" "TMe-755"  "TMe-2213"
#> [25] "TMe-803"  "TMe-863"  "TMe-1004" "TMe-1500" "TMe-667"  "TMe-1003"
#> [31] "TMe-712"  "TMe-2124" "TMe-2753" "TMe-588"  "TMe-603"  "TMe-2750"
#> [37] "TMe-2290" "TMe-348"  "TMe-892"  "TMe-3329"
#> 
#> [[6]]
#>  [1] "TMe-1392" "TMe-1302" "TMe-1217" "TMe-1076" "TMe-1608" "TMe-854" 
#>  [7] "TMe-1413" "TMe-1481" "TMe-751"  "TMe-1124" "TMe-2791" "TMe-620" 
#> [13] "TMe-968"  "TMe-693"  "TMe-2543" "TMe-1661" "TMe-3177"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(randomsel_sum_margalef,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-3694" "TMe-3236" "TMe-1086"
#>  [7] "TMe-865"  "TMe-2967" "TMe-3481" "TMe-3051" "TMe-2975" "TMe-1922"
#> [13] "TMe-3252" "TMe-3262" "TMe-937"  "TMe-1218" "TMe-2810" "TMe-3112"
#> [19] "TMe-2027" "TMe-778"  "TMe-1717" "TMe-3726" "TMe-1581" "TMe-2066"
#> [25] "TMe-469"  "TMe-500"  "TMe-952"  "TMe-1170" "TMe-1532" "TMe-1914"
#> [31] "TMe-1451"
#> 
#> [[2]]
#>  [1] "TMe-2757" "TMe-3530" "TMe-1754" "TMe-2329" "TMe-1385" "TMe-2352"
#>  [7] "TMe-339"  "TMe-74"   "TMe-369"  "TMe-3547" "TMe-3366" "TMe-2033"
#> [13] "TMe-40"   "TMe-3200" "TMe-2412" "TMe-2715" "TMe-3800" "TMe-681" 
#> [19] "TMe-2951" "TMe-196"  "TMe-2995" "TMe-674"  "TMe-1831" "TMe-2997"
#> [25] "TMe-796"  "TMe-2568" "TMe-3766" "TMe-1505" "TMe-1668" "TMe-1732"
#> [31] "TMe-3495"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-3133" "TMe-3397" "TMe-234"  "TMe-1804" "TMe-2161" "TMe-1738"
#> [13] "TMe-3100" "TMe-3336" "TMe-1787" "TMe-1939" "TMe-3644" "TMe-123" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-2458" "TMe-3255" "TMe-3327"
#>  [7] "TMe-1376" "TMe-698"  "TMe-1434" "TMe-280"  "TMe-3428" "TMe-3527"
#> [13] "TMe-241"  "TMe-737"  "TMe-1020" "TMe-259"  "TMe-279"  "TMe-25"  
#> [19] "TMe-619"  "TMe-3273" "TMe-1988" "TMe-1148" "TMe-2956" "TMe-1297"
#> [25] "TMe-699"  "TMe-1987" "TMe-1155" "TMe-2375" "TMe-191"  "TMe-1580"
#> [31] "TMe-1144" "TMe-3072" "TMe-57"   "TMe-2318" "TMe-1016"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-723"  "TMe-419" 
#>  [7] "TMe-3329" "TMe-1004" "TMe-2853" "TMe-1003" "TMe-423"  "TMe-1427"
#> [13] "TMe-588"  "TMe-574"  "TMe-629"  "TMe-688"  "TMe-333"  "TMe-2124"
#> [19] "TMe-1375" "TMe-769"  "TMe-344"  "TMe-1220" "TMe-532"  "TMe-745" 
#> [25] "TMe-1534" "TMe-2151" "TMe-1293" "TMe-332"  "TMe-247"  "TMe-1012"
#> [31] "TMe-1290" "TMe-863"  "TMe-439"  "TMe-682"  "TMe-1131" "TMe-487" 
#> [37] "TMe-997"  "TMe-1880" "TMe-1268" "TMe-1760"
#> 
#> [[6]]
#>  [1] "TMe-625"  "TMe-1124" "TMe-2791" "TMe-1566" "TMe-2543" "TMe-1392"
#>  [7] "TMe-1816" "TMe-693"  "TMe-751"  "TMe-1608" "TMe-1900" "TMe-1302"
#> [13] "TMe-1756" "TMe-1592" "TMe-1481" "TMe-1239" "TMe-1062"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_mean_richness,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-3694" "TMe-3236" "TMe-1086"
#>  [7] "TMe-865"  "TMe-2967" "TMe-3481" "TMe-3051" "TMe-2975" "TMe-1922"
#> [13] "TMe-3252" "TMe-3262" "TMe-937"  "TMe-1218" "TMe-2810" "TMe-3112"
#> [19] "TMe-2027" "TMe-778"  "TMe-1717" "TMe-3726" "TMe-1581" "TMe-2066"
#> [25] "TMe-469"  "TMe-500"  "TMe-952"  "TMe-1170" "TMe-1532" "TMe-1914"
#> [31] "TMe-1451"
#> 
#> [[2]]
#>  [1] "TMe-40"   "TMe-2352" "TMe-2329" "TMe-3530" "TMe-2995" "TMe-681" 
#>  [7] "TMe-2568" "TMe-960"  "TMe-369"  "TMe-74"   "TMe-3366" "TMe-2033"
#> [13] "TMe-3547" "TMe-3200" "TMe-2412" "TMe-2715" "TMe-3800" "TMe-2951"
#> [19] "TMe-196"  "TMe-339"  "TMe-674"  "TMe-1831" "TMe-2997" "TMe-796" 
#> [25] "TMe-3766" "TMe-1505" "TMe-1668" "TMe-1732" "TMe-3495" "TMe-2757"
#> [31] "TMe-3557"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-3133" "TMe-3397" "TMe-234"  "TMe-1804" "TMe-2161" "TMe-1738"
#> [13] "TMe-3100" "TMe-3336" "TMe-1787" "TMe-1939" "TMe-3644" "TMe-123" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-2807" "TMe-3255" "TMe-698" 
#>  [7] "TMe-266"  "TMe-1434" "TMe-3248" "TMe-1020" "TMe-3527" "TMe-1579"
#> [13] "TMe-241"  "TMe-737"  "TMe-279"  "TMe-25"   "TMe-3428" "TMe-619" 
#> [19] "TMe-3273" "TMe-280"  "TMe-1148" "TMe-2956" "TMe-1297" "TMe-699" 
#> [25] "TMe-259"  "TMe-1987" "TMe-1155" "TMe-2375" "TMe-191"  "TMe-1580"
#> [31] "TMe-1144" "TMe-3072" "TMe-1376" "TMe-57"   "TMe-2318"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-723"  "TMe-419" 
#>  [7] "TMe-3329" "TMe-1004" "TMe-2853" "TMe-1003" "TMe-423"  "TMe-1427"
#> [13] "TMe-588"  "TMe-574"  "TMe-629"  "TMe-688"  "TMe-333"  "TMe-2124"
#> [19] "TMe-1375" "TMe-769"  "TMe-344"  "TMe-1220" "TMe-532"  "TMe-745" 
#> [25] "TMe-1534" "TMe-2151" "TMe-1293" "TMe-332"  "TMe-247"  "TMe-1012"
#> [31] "TMe-1290" "TMe-863"  "TMe-439"  "TMe-682"  "TMe-1131" "TMe-487" 
#> [37] "TMe-997"  "TMe-1880" "TMe-1268" "TMe-1760"
#> 
#> [[6]]
#>  [1] "TMe-1481" "TMe-1403" "TMe-693"  "TMe-1392" "TMe-1646" "TMe-1816"
#>  [7] "TMe-2791" "TMe-1566" "TMe-1302" "TMe-1756" "TMe-1124" "TMe-751" 
#> [13] "TMe-1900" "TMe-1592" "TMe-1608" "TMe-1239" "TMe-1062"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_sum_richness,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2453" "TMe-2785"
#>  [7] "TMe-1096" "TMe-778"  "TMe-3112" "TMe-41"   "TMe-3130" "TMe-1086"
#> [13] "TMe-2027" "TMe-1190" "TMe-3481" "TMe-3115" "TMe-3262" "TMe-1922"
#> [19] "TMe-1218" "TMe-3252" "TMe-3437" "TMe-1869" "TMe-2237" "TMe-566" 
#> [25] "TMe-2975" "TMe-3398" "TMe-3111" "TMe-2810" "TMe-569"  "TMe-937" 
#> [31] "TMe-3236"
#> 
#> [[2]]
#>  [1] "TMe-796"  "TMe-3805" "TMe-171"  "TMe-3200" "TMe-3258" "TMe-369" 
#>  [7] "TMe-40"   "TMe-2903" "TMe-1385" "TMe-2329" "TMe-3766" "TMe-2952"
#> [13] "TMe-2033" "TMe-1754" "TMe-2568" "TMe-539"  "TMe-74"   "TMe-3547"
#> [19] "TMe-674"  "TMe-1668" "TMe-3366" "TMe-3447" "TMe-3101" "TMe-2021"
#> [25] "TMe-3530" "TMe-2352" "TMe-2211" "TMe-2000" "TMe-2257" "TMe-3284"
#> [31] "TMe-3093"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-234"  "TMe-3397" "TMe-3569" "TMe-1804" "TMe-3100" "TMe-425" 
#> [13] "TMe-3631" "TMe-2161" "TMe-1230" "TMe-1819" "TMe-261"  "TMe-3148"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-3248" "TMe-698" 
#>  [7] "TMe-1434" "TMe-3527" "TMe-956"  "TMe-241"  "TMe-3542" "TMe-1020"
#> [13] "TMe-2958" "TMe-1988" "TMe-2924" "TMe-3273" "TMe-619"  "TMe-18"  
#> [19] "TMe-2971" "TMe-3428" "TMe-812"  "TMe-1377" "TMe-280"  "TMe-266" 
#> [25] "TMe-3189" "TMe-3242" "TMe-2956" "TMe-1700" "TMe-1776" "TMe-1123"
#> [31] "TMe-25"   "TMe-1376" "TMe-279"  "TMe-427"  "TMe-3390"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-419"  "TMe-723" 
#>  [7] "TMe-3329" "TMe-2003" "TMe-2853" "TMe-877"  "TMe-2290" "TMe-423" 
#> [13] "TMe-1375" "TMe-2355" "TMe-627"  "TMe-755"  "TMe-532"  "TMe-1220"
#> [19] "TMe-2907" "TMe-1541" "TMe-167"  "TMe-1004" "TMe-2425" "TMe-1500"
#> [25] "TMe-816"  "TMe-1427" "TMe-2750" "TMe-730"  "TMe-667"  "TMe-487" 
#> [31] "TMe-1534" "TMe-344"  "TMe-98"   "TMe-1196" "TMe-1934" "TMe-1290"
#> [37] "TMe-823"  "TMe-2855" "TMe-2192" "TMe-439" 
#> 
#> [[6]]
#>  [1] "TMe-661"  "TMe-310"  "TMe-1035" "TMe-1403" "TMe-1392" "TMe-2791"
#>  [7] "TMe-1124" "TMe-1566" "TMe-1816" "TMe-693"  "TMe-2983" "TMe-2818"
#> [13] "TMe-1661" "TMe-751"  "TMe-2543" "TMe-598"  "TMe-1900"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_mean_shannon,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2453" "TMe-2785"
#>  [7] "TMe-1096" "TMe-778"  "TMe-3112" "TMe-41"   "TMe-3130" "TMe-1086"
#> [13] "TMe-2027" "TMe-1190" "TMe-3481" "TMe-3115" "TMe-3262" "TMe-1922"
#> [19] "TMe-1218" "TMe-3252" "TMe-3437" "TMe-1869" "TMe-2237" "TMe-566" 
#> [25] "TMe-2975" "TMe-3398" "TMe-3111" "TMe-2810" "TMe-569"  "TMe-937" 
#> [31] "TMe-3236"
#> 
#> [[2]]
#>  [1] "TMe-251"  "TMe-2000" "TMe-3258" "TMe-2568" "TMe-369"  "TMe-2033"
#>  [7] "TMe-2211" "TMe-2329" "TMe-3200" "TMe-40"   "TMe-3805" "TMe-1385"
#> [13] "TMe-1754" "TMe-3766" "TMe-2952" "TMe-3366" "TMe-3530" "TMe-171" 
#> [19] "TMe-3547" "TMe-3101" "TMe-2352" "TMe-2021" "TMe-674"  "TMe-539" 
#> [25] "TMe-74"   "TMe-1668" "TMe-3447" "TMe-2903" "TMe-3093" "TMe-3284"
#> [31] "TMe-2257"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-234"  "TMe-3397" "TMe-3569" "TMe-1804" "TMe-3100" "TMe-425" 
#> [13] "TMe-3631" "TMe-2161" "TMe-1230" "TMe-1819" "TMe-261"  "TMe-3148"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-3248" "TMe-698" 
#>  [7] "TMe-1434" "TMe-3527" "TMe-956"  "TMe-241"  "TMe-3542" "TMe-1020"
#> [13] "TMe-2958" "TMe-1988" "TMe-2924" "TMe-3273" "TMe-619"  "TMe-18"  
#> [19] "TMe-2971" "TMe-3428" "TMe-812"  "TMe-1377" "TMe-280"  "TMe-266" 
#> [25] "TMe-3189" "TMe-3242" "TMe-2956" "TMe-1700" "TMe-1776" "TMe-1123"
#> [31] "TMe-25"   "TMe-1376" "TMe-279"  "TMe-427"  "TMe-3390"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-419"  "TMe-723" 
#>  [7] "TMe-3329" "TMe-2003" "TMe-2853" "TMe-877"  "TMe-2290" "TMe-423" 
#> [13] "TMe-1375" "TMe-2355" "TMe-627"  "TMe-755"  "TMe-532"  "TMe-1220"
#> [19] "TMe-2907" "TMe-1541" "TMe-167"  "TMe-1004" "TMe-2425" "TMe-1500"
#> [25] "TMe-816"  "TMe-1427" "TMe-2750" "TMe-730"  "TMe-667"  "TMe-487" 
#> [31] "TMe-1534" "TMe-344"  "TMe-98"   "TMe-1196" "TMe-1934" "TMe-1290"
#> [37] "TMe-823"  "TMe-2855" "TMe-2192" "TMe-439" 
#> 
#> [[6]]
#>  [1] "TMe-661"  "TMe-310"  "TMe-1035" "TMe-1403" "TMe-1392" "TMe-2791"
#>  [7] "TMe-1124" "TMe-1566" "TMe-1816" "TMe-693"  "TMe-2983" "TMe-2818"
#> [13] "TMe-1661" "TMe-751"  "TMe-2543" "TMe-598"  "TMe-1900"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_sum_shannon,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-566" 
#>  [7] "TMe-2453" "TMe-41"   "TMe-3252" "TMe-1096" "TMe-500"  "TMe-3112"
#> [13] "TMe-1190" "TMe-3437" "TMe-2237" "TMe-3115" "TMe-1922" "TMe-778" 
#> [19] "TMe-1086" "TMe-3685" "TMe-2027" "TMe-2810" "TMe-1960" "TMe-1589"
#> [25] "TMe-3262" "TMe-3398" "TMe-1869" "TMe-569"  "TMe-3236" "TMe-937" 
#> [31] "TMe-2975"
#> 
#> [[2]]
#>  [1] "TMe-3258" "TMe-369"  "TMe-3805" "TMe-1754" "TMe-3366" "TMe-2033"
#>  [7] "TMe-2952" "TMe-674"  "TMe-2329" "TMe-3200" "TMe-40"   "TMe-1385"
#> [13] "TMe-3766" "TMe-3101" "TMe-171"  "TMe-3447" "TMe-2568" "TMe-3530"
#> [19] "TMe-2352" "TMe-539"  "TMe-2021" "TMe-3284" "TMe-2257" "TMe-3547"
#> [25] "TMe-2000" "TMe-2211" "TMe-1732" "TMe-251"  "TMe-2903" "TMe-796" 
#> [31] "TMe-3093"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1863" "TMe-3148" "TMe-2161" "TMe-3631" "TMe-425"  "TMe-1804"
#> [13] "TMe-1819" "TMe-261"  "TMe-3100" "TMe-2374" "TMe-1230" "TMe-3397"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-241"  "TMe-2924" "TMe-812" 
#>  [7] "TMe-3542" "TMe-1377" "TMe-2971" "TMe-3273" "TMe-619"  "TMe-3527"
#> [13] "TMe-1434" "TMe-698"  "TMe-3242" "TMe-280"  "TMe-1020" "TMe-2958"
#> [19] "TMe-3189" "TMe-1776" "TMe-1988" "TMe-3390" "TMe-3255" "TMe-1376"
#> [25] "TMe-1297" "TMe-3198" "TMe-1350" "TMe-3428" "TMe-2956" "TMe-1144"
#> [31] "TMe-25"   "TMe-372"  "TMe-259"  "TMe-266"  "TMe-1579"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-2290"
#>  [7] "TMe-2355" "TMe-2003" "TMe-755"  "TMe-419"  "TMe-3329" "TMe-627" 
#> [13] "TMe-1375" "TMe-2853" "TMe-423"  "TMe-1541" "TMe-2750" "TMe-344" 
#> [19] "TMe-823"  "TMe-1534" "TMe-730"  "TMe-2192" "TMe-1500" "TMe-532" 
#> [25] "TMe-487"  "TMe-2907" "TMe-877"  "TMe-1196" "TMe-167"  "TMe-1159"
#> [31] "TMe-795"  "TMe-98"   "TMe-667"  "TMe-2855" "TMe-816"  "TMe-224" 
#> [37] "TMe-1934" "TMe-1265" "TMe-1427" "TMe-688" 
#> 
#> [[6]]
#>  [1] "TMe-1445" "TMe-1403" "TMe-2791" "TMe-1566" "TMe-2983" "TMe-1816"
#>  [7] "TMe-2818" "TMe-693"  "TMe-1124" "TMe-1661" "TMe-1035" "TMe-2543"
#> [13] "TMe-751"  "TMe-598"  "TMe-1392" "TMe-310"  "TMe-1992"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_mean_simpson,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-1190"
#>  [7] "TMe-2453" "TMe-41"   "TMe-1589" "TMe-3262" "TMe-1096" "TMe-3481"
#> [13] "TMe-1086" "TMe-3115" "TMe-3437" "TMe-778"  "TMe-2027" "TMe-3252"
#> [19] "TMe-566"  "TMe-3112" "TMe-1922" "TMe-1869" "TMe-2810" "TMe-3685"
#> [25] "TMe-2237" "TMe-1218" "TMe-569"  "TMe-3111" "TMe-2975" "TMe-300" 
#> [31] "TMe-3398"
#> 
#> [[2]]
#>  [1] "TMe-1754" "TMe-3366" "TMe-369"  "TMe-2352" "TMe-2033" "TMe-2952"
#>  [7] "TMe-3766" "TMe-2329" "TMe-40"   "TMe-1385" "TMe-2568" "TMe-2211"
#> [13] "TMe-3530" "TMe-171"  "TMe-3200" "TMe-539"  "TMe-3101" "TMe-3805"
#> [19] "TMe-674"  "TMe-2257" "TMe-2000" "TMe-3284" "TMe-2021" "TMe-3447"
#> [25] "TMe-3547" "TMe-251"  "TMe-2903" "TMe-1831" "TMe-3093" "TMe-1732"
#> [31] "TMe-3258"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1863" "TMe-3148" "TMe-2161" "TMe-3631" "TMe-425"  "TMe-1804"
#> [13] "TMe-2374" "TMe-3100" "TMe-261"  "TMe-1230" "TMe-1819" "TMe-3397"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-241"  "TMe-2971" "TMe-1377"
#>  [7] "TMe-2924" "TMe-3255" "TMe-698"  "TMe-1434" "TMe-3527" "TMe-1020"
#> [13] "TMe-3542" "TMe-3273" "TMe-3242" "TMe-619"  "TMe-3390" "TMe-280" 
#> [19] "TMe-812"  "TMe-3189" "TMe-2958" "TMe-1988" "TMe-1376" "TMe-1776"
#> [25] "TMe-1297" "TMe-3198" "TMe-1350" "TMe-3428" "TMe-2956" "TMe-1144"
#> [31] "TMe-25"   "TMe-372"  "TMe-259"  "TMe-266"  "TMe-1579"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-2290"
#>  [7] "TMe-2355" "TMe-2003" "TMe-755"  "TMe-419"  "TMe-3329" "TMe-627" 
#> [13] "TMe-1375" "TMe-2853" "TMe-423"  "TMe-1541" "TMe-2750" "TMe-344" 
#> [19] "TMe-823"  "TMe-1534" "TMe-730"  "TMe-2192" "TMe-1500" "TMe-532" 
#> [25] "TMe-487"  "TMe-2907" "TMe-877"  "TMe-1196" "TMe-167"  "TMe-1159"
#> [31] "TMe-795"  "TMe-98"   "TMe-667"  "TMe-2855" "TMe-816"  "TMe-224" 
#> [37] "TMe-1934" "TMe-1265" "TMe-1427" "TMe-688" 
#> 
#> [[6]]
#>  [1] "TMe-1403" "TMe-1445" "TMe-1816" "TMe-1566" "TMe-2983" "TMe-2791"
#>  [7] "TMe-693"  "TMe-2543" "TMe-1124" "TMe-1035" "TMe-1661" "TMe-598" 
#> [13] "TMe-751"  "TMe-2818" "TMe-310"  "TMe-1392" "TMe-1992"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_sum_simpson,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-2453"
#>  [7] "TMe-1190" "TMe-41"   "TMe-3685" "TMe-2027" "TMe-3719" "TMe-1096"
#> [13] "TMe-778"  "TMe-1086" "TMe-3115" "TMe-3262" "TMe-566"  "TMe-3112"
#> [19] "TMe-1922" "TMe-3437" "TMe-3539" "TMe-3236" "TMe-569"  "TMe-2237"
#> [25] "TMe-2810" "TMe-300"  "TMe-1589" "TMe-3398" "TMe-1869" "TMe-865" 
#> [31] "TMe-2975"
#> 
#> [[2]]
#>  [1] "TMe-2611" "TMe-369"  "TMe-3805" "TMe-171"  "TMe-2952" "TMe-2033"
#>  [7] "TMe-3200" "TMe-2329" "TMe-1385" "TMe-1754" "TMe-2568" "TMe-3766"
#> [13] "TMe-40"   "TMe-1831" "TMe-3530" "TMe-3366" "TMe-2211" "TMe-3101"
#> [19] "TMe-2352" "TMe-539"  "TMe-2021" "TMe-674"  "TMe-3284" "TMe-2257"
#> [25] "TMe-3258" "TMe-2000" "TMe-2903" "TMe-251"  "TMe-3447" "TMe-3093"
#> [31] "TMe-796" 
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1230" "TMe-1863" "TMe-1819" "TMe-261"  "TMe-3148" "TMe-1804"
#> [13] "TMe-425"  "TMe-3631" "TMe-2161" "TMe-2374" "TMe-3397" "TMe-3100"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-1144" "TMe-1434"
#>  [7] "TMe-3527" "TMe-698"  "TMe-2971" "TMe-241"  "TMe-3409" "TMe-1988"
#> [13] "TMe-2924" "TMe-3273" "TMe-1376" "TMe-2958" "TMe-1020" "TMe-3542"
#> [19] "TMe-3242" "TMe-280"  "TMe-812"  "TMe-1776" "TMe-3189" "TMe-2956"
#> [25] "TMe-3198" "TMe-1377" "TMe-3428" "TMe-353"  "TMe-3390" "TMe-1297"
#> [31] "TMe-1579" "TMe-372"  "TMe-25"   "TMe-1350" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-2290"
#>  [7] "TMe-2355" "TMe-2003" "TMe-755"  "TMe-419"  "TMe-3329" "TMe-627" 
#> [13] "TMe-1220" "TMe-2192" "TMe-2853" "TMe-423"  "TMe-1375" "TMe-532" 
#> [19] "TMe-877"  "TMe-1500" "TMe-2750" "TMe-730"  "TMe-167"  "TMe-2907"
#> [25] "TMe-816"  "TMe-667"  "TMe-1541" "TMe-487"  "TMe-1534" "TMe-344" 
#> [31] "TMe-1196" "TMe-920"  "TMe-98"   "TMe-1401" "TMe-1159" "TMe-823" 
#> [37] "TMe-1934" "TMe-224"  "TMe-1004" "TMe-1295"
#> 
#> [[6]]
#>  [1] "TMe-2791" "TMe-598"  "TMe-1403" "TMe-1566" "TMe-1816" "TMe-1992"
#>  [7] "TMe-693"  "TMe-2543" "TMe-1035" "TMe-751"  "TMe-2818" "TMe-1661"
#> [13] "TMe-1124" "TMe-1392" "TMe-2983" "TMe-310"  "TMe-1503"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_mean_mcintosh,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-2453"
#>  [7] "TMe-1190" "TMe-41"   "TMe-3685" "TMe-2027" "TMe-3719" "TMe-1096"
#> [13] "TMe-778"  "TMe-1086" "TMe-3115" "TMe-3262" "TMe-566"  "TMe-3112"
#> [19] "TMe-1922" "TMe-3437" "TMe-3539" "TMe-3236" "TMe-569"  "TMe-2237"
#> [25] "TMe-2810" "TMe-300"  "TMe-1589" "TMe-3398" "TMe-1869" "TMe-865" 
#> [31] "TMe-2975"
#> 
#> [[2]]
#>  [1] "TMe-2211" "TMe-2329" "TMe-3530" "TMe-2352" "TMe-40"   "TMe-3766"
#>  [7] "TMe-369"  "TMe-1385" "TMe-2952" "TMe-2033" "TMe-1754" "TMe-3200"
#> [13] "TMe-171"  "TMe-2568" "TMe-3101" "TMe-3547" "TMe-1668" "TMe-3805"
#> [19] "TMe-674"  "TMe-3366" "TMe-3447" "TMe-539"  "TMe-3284" "TMe-2021"
#> [25] "TMe-251"  "TMe-2257" "TMe-2903" "TMe-2000" "TMe-3093" "TMe-796" 
#> [31] "TMe-3258"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1230" "TMe-1863" "TMe-1819" "TMe-261"  "TMe-3148" "TMe-1804"
#> [13] "TMe-425"  "TMe-3631" "TMe-2161" "TMe-2374" "TMe-3397" "TMe-3100"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-1144" "TMe-1434"
#>  [7] "TMe-3527" "TMe-698"  "TMe-2971" "TMe-241"  "TMe-3409" "TMe-1988"
#> [13] "TMe-2924" "TMe-3273" "TMe-1376" "TMe-2958" "TMe-1020" "TMe-3542"
#> [19] "TMe-3242" "TMe-280"  "TMe-812"  "TMe-1776" "TMe-3189" "TMe-2956"
#> [25] "TMe-3198" "TMe-1377" "TMe-3428" "TMe-353"  "TMe-3390" "TMe-1297"
#> [31] "TMe-1579" "TMe-372"  "TMe-25"   "TMe-1350" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-2290"
#>  [7] "TMe-2355" "TMe-2003" "TMe-755"  "TMe-419"  "TMe-3329" "TMe-627" 
#> [13] "TMe-1220" "TMe-2192" "TMe-2853" "TMe-423"  "TMe-1375" "TMe-532" 
#> [19] "TMe-877"  "TMe-1500" "TMe-2750" "TMe-730"  "TMe-167"  "TMe-2907"
#> [25] "TMe-816"  "TMe-667"  "TMe-1541" "TMe-487"  "TMe-1534" "TMe-344" 
#> [31] "TMe-1196" "TMe-920"  "TMe-98"   "TMe-1401" "TMe-1159" "TMe-823" 
#> [37] "TMe-1934" "TMe-224"  "TMe-1004" "TMe-1295"
#> 
#> [[6]]
#>  [1] "TMe-598"  "TMe-1124" "TMe-2791" "TMe-1566" "TMe-1816" "TMe-1403"
#>  [7] "TMe-1992" "TMe-2983" "TMe-1661" "TMe-1035" "TMe-2818" "TMe-751" 
#> [13] "TMe-2543" "TMe-693"  "TMe-1392" "TMe-310"  "TMe-1503"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_sum_mcintosh,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-2453"
#>  [7] "TMe-1190" "TMe-41"   "TMe-3481" "TMe-2027" "TMe-3130" "TMe-1096"
#> [13] "TMe-778"  "TMe-1086" "TMe-3115" "TMe-3112" "TMe-1922" "TMe-3437"
#> [19] "TMe-566"  "TMe-3262" "TMe-2237" "TMe-3539" "TMe-2810" "TMe-1869"
#> [25] "TMe-3398" "TMe-3252" "TMe-2975" "TMe-3111" "TMe-569"  "TMe-1218"
#> [31] "TMe-300" 
#> 
#> [[2]]
#>  [1] "TMe-3258" "TMe-2033" "TMe-369"  "TMe-3766" "TMe-2568" "TMe-1754"
#>  [7] "TMe-2952" "TMe-1385" "TMe-2903" "TMe-2329" "TMe-3200" "TMe-40"  
#> [13] "TMe-171"  "TMe-3530" "TMe-3366" "TMe-3805" "TMe-3101" "TMe-2352"
#> [19] "TMe-3547" "TMe-1668" "TMe-539"  "TMe-674"  "TMe-251"  "TMe-2021"
#> [25] "TMe-74"   "TMe-3447" "TMe-3093" "TMe-2000" "TMe-2257" "TMe-796" 
#> [31] "TMe-2211"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1863" "TMe-3397" "TMe-3569" "TMe-3100" "TMe-1804" "TMe-425" 
#> [13] "TMe-3631" "TMe-1819" "TMe-1230" "TMe-2161" "TMe-261"  "TMe-3148"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-1434" "TMe-1700"
#>  [7] "TMe-698"  "TMe-3527" "TMe-2971" "TMe-241"  "TMe-1155" "TMe-2924"
#> [13] "TMe-1988" "TMe-3273" "TMe-1020" "TMe-2958" "TMe-18"   "TMe-372" 
#> [19] "TMe-280"  "TMe-812"  "TMe-3542" "TMe-1776" "TMe-3428" "TMe-1376"
#> [25] "TMe-956"  "TMe-3189" "TMe-3242" "TMe-1297" "TMe-2956" "TMe-266" 
#> [31] "TMe-25"   "TMe-1123" "TMe-3390" "TMe-1579" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-3329"
#>  [7] "TMe-419"  "TMe-2853" "TMe-2003" "TMe-755"  "TMe-423"  "TMe-2290"
#> [13] "TMe-1375" "TMe-2355" "TMe-627"  "TMe-1220" "TMe-532"  "TMe-1500"
#> [19] "TMe-877"  "TMe-2750" "TMe-2907" "TMe-167"  "TMe-816"  "TMe-1004"
#> [25] "TMe-1427" "TMe-487"  "TMe-667"  "TMe-730"  "TMe-1534" "TMe-1196"
#> [31] "TMe-344"  "TMe-98"   "TMe-2855" "TMe-823"  "TMe-688"  "TMe-2192"
#> [37] "TMe-362"  "TMe-1290" "TMe-1541" "TMe-1934"
#> 
#> [[6]]
#>  [1] "TMe-1392" "TMe-693"  "TMe-1403" "TMe-1816" "TMe-1124" "TMe-2791"
#>  [7] "TMe-1566" "TMe-2983" "TMe-1661" "TMe-2818" "TMe-1035" "TMe-751" 
#> [13] "TMe-2543" "TMe-310"  "TMe-661"  "TMe-598"  "TMe-1900"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_mean_brillouin,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-2453"
#>  [7] "TMe-1190" "TMe-41"   "TMe-3481" "TMe-2027" "TMe-3130" "TMe-1096"
#> [13] "TMe-778"  "TMe-1086" "TMe-3115" "TMe-3112" "TMe-1922" "TMe-3437"
#> [19] "TMe-566"  "TMe-3262" "TMe-2237" "TMe-3539" "TMe-2810" "TMe-1869"
#> [25] "TMe-3398" "TMe-3252" "TMe-2975" "TMe-3111" "TMe-569"  "TMe-1218"
#> [31] "TMe-300" 
#> 
#> [[2]]
#>  [1] "TMe-674"  "TMe-1754" "TMe-2033" "TMe-2952" "TMe-369"  "TMe-3366"
#>  [7] "TMe-2211" "TMe-2329" "TMe-3200" "TMe-3805" "TMe-40"   "TMe-1385"
#> [13] "TMe-3766" "TMe-3101" "TMe-3530" "TMe-2021" "TMe-3547" "TMe-171" 
#> [19] "TMe-2568" "TMe-2352" "TMe-539"  "TMe-74"   "TMe-3447" "TMe-1668"
#> [25] "TMe-2000" "TMe-2257" "TMe-3093" "TMe-251"  "TMe-2903" "TMe-796" 
#> [31] "TMe-3258"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1863" "TMe-3397" "TMe-3569" "TMe-3100" "TMe-1804" "TMe-425" 
#> [13] "TMe-3631" "TMe-1819" "TMe-1230" "TMe-2161" "TMe-261"  "TMe-3148"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-1434" "TMe-1700"
#>  [7] "TMe-698"  "TMe-3527" "TMe-2971" "TMe-241"  "TMe-1155" "TMe-2924"
#> [13] "TMe-1988" "TMe-3273" "TMe-1020" "TMe-2958" "TMe-18"   "TMe-372" 
#> [19] "TMe-280"  "TMe-812"  "TMe-3542" "TMe-1776" "TMe-3428" "TMe-1376"
#> [25] "TMe-956"  "TMe-3189" "TMe-3242" "TMe-1297" "TMe-2956" "TMe-266" 
#> [31] "TMe-25"   "TMe-1123" "TMe-3390" "TMe-1579" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-3329"
#>  [7] "TMe-419"  "TMe-2853" "TMe-2003" "TMe-755"  "TMe-423"  "TMe-2290"
#> [13] "TMe-1375" "TMe-2355" "TMe-627"  "TMe-1220" "TMe-532"  "TMe-1500"
#> [19] "TMe-877"  "TMe-2750" "TMe-2907" "TMe-167"  "TMe-816"  "TMe-1004"
#> [25] "TMe-1427" "TMe-487"  "TMe-667"  "TMe-730"  "TMe-1534" "TMe-1196"
#> [31] "TMe-344"  "TMe-98"   "TMe-2855" "TMe-823"  "TMe-688"  "TMe-2192"
#> [37] "TMe-362"  "TMe-1290" "TMe-1541" "TMe-1934"
#> 
#> [[6]]
#>  [1] "TMe-1445" "TMe-1992" "TMe-1661" "TMe-2983" "TMe-1566" "TMe-2791"
#>  [7] "TMe-2543" "TMe-693"  "TMe-1816" "TMe-1124" "TMe-1392" "TMe-598" 
#> [13] "TMe-1403" "TMe-751"  "TMe-310"  "TMe-1035" "TMe-2818"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_sum_brillouin,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-3694" "TMe-3236" "TMe-1086"
#>  [7] "TMe-865"  "TMe-41"   "TMe-3051" "TMe-2785" "TMe-2975" "TMe-1922"
#> [13] "TMe-3252" "TMe-3262" "TMe-937"  "TMe-1218" "TMe-3112" "TMe-2027"
#> [19] "TMe-778"  "TMe-1717" "TMe-3726" "TMe-1581" "TMe-2066" "TMe-469" 
#> [25] "TMe-500"  "TMe-952"  "TMe-1170" "TMe-1532" "TMe-1914" "TMe-1451"
#> [31] "TMe-1117"
#> 
#> [[2]]
#>  [1] "TMe-960"  "TMe-3766" "TMe-2952" "TMe-2033" "TMe-796"  "TMe-40"  
#>  [7] "TMe-369"  "TMe-1385" "TMe-3800" "TMe-3547" "TMe-74"   "TMe-3366"
#> [13] "TMe-3200" "TMe-2412" "TMe-2715" "TMe-681"  "TMe-2951" "TMe-196" 
#> [19] "TMe-339"  "TMe-2995" "TMe-674"  "TMe-2352" "TMe-1831" "TMe-2997"
#> [25] "TMe-2568" "TMe-1505" "TMe-2329" "TMe-1668" "TMe-1732" "TMe-3495"
#> [31] "TMe-2757"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-3133" "TMe-3397" "TMe-234"  "TMe-1804" "TMe-2161" "TMe-1738"
#> [13] "TMe-3336" "TMe-3100" "TMe-1787" "TMe-1939" "TMe-3644" "TMe-123" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-2458" "TMe-3255" "TMe-3327"
#>  [7] "TMe-1376" "TMe-698"  "TMe-1434" "TMe-280"  "TMe-3428" "TMe-3527"
#> [13] "TMe-241"  "TMe-737"  "TMe-1020" "TMe-3273" "TMe-259"  "TMe-279" 
#> [19] "TMe-25"   "TMe-619"  "TMe-1988" "TMe-1148" "TMe-2956" "TMe-1297"
#> [25] "TMe-699"  "TMe-1987" "TMe-1155" "TMe-2375" "TMe-191"  "TMe-1580"
#> [31] "TMe-1144" "TMe-3072" "TMe-57"   "TMe-2318" "TMe-1016"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-723"  "TMe-419" 
#>  [7] "TMe-3329" "TMe-1004" "TMe-2853" "TMe-1003" "TMe-423"  "TMe-1427"
#> [13] "TMe-588"  "TMe-574"  "TMe-629"  "TMe-688"  "TMe-333"  "TMe-2124"
#> [19] "TMe-1375" "TMe-769"  "TMe-344"  "TMe-1220" "TMe-532"  "TMe-745" 
#> [25] "TMe-1534" "TMe-2151" "TMe-1293" "TMe-332"  "TMe-247"  "TMe-1012"
#> [31] "TMe-1290" "TMe-863"  "TMe-439"  "TMe-682"  "TMe-1131" "TMe-487" 
#> [37] "TMe-997"  "TMe-1880" "TMe-1268" "TMe-1760"
#> 
#> [[6]]
#>  [1] "TMe-809"  "TMe-1403" "TMe-2791" "TMe-1062" "TMe-1816" "TMe-693" 
#>  [7] "TMe-2983" "TMe-1302" "TMe-1392" "TMe-1756" "TMe-1124" "TMe-751" 
#> [13] "TMe-1900" "TMe-1592" "TMe-1566" "TMe-1481" "TMe-1608"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_mean_margalef,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-3694" "TMe-3236" "TMe-1086"
#>  [7] "TMe-3481" "TMe-865"  "TMe-2967" "TMe-3051" "TMe-2975" "TMe-1922"
#> [13] "TMe-3252" "TMe-3262" "TMe-2810" "TMe-937"  "TMe-1218" "TMe-2027"
#> [19] "TMe-778"  "TMe-3112" "TMe-1717" "TMe-3726" "TMe-1581" "TMe-2066"
#> [25] "TMe-469"  "TMe-500"  "TMe-952"  "TMe-1170" "TMe-1532" "TMe-1914"
#> [31] "TMe-1451"
#> 
#> [[2]]
#>  [1] "TMe-2997" "TMe-3805" "TMe-171"  "TMe-3200" "TMe-1831" "TMe-40"  
#>  [7] "TMe-369"  "TMe-674"  "TMe-539"  "TMe-74"   "TMe-2329" "TMe-2033"
#> [13] "TMe-3547" "TMe-2412" "TMe-2715" "TMe-3800" "TMe-681"  "TMe-2951"
#> [19] "TMe-196"  "TMe-339"  "TMe-2995" "TMe-2352" "TMe-796"  "TMe-2568"
#> [25] "TMe-3766" "TMe-1505" "TMe-1668" "TMe-1732" "TMe-3495" "TMe-960" 
#> [31] "TMe-2757"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-3133" "TMe-3397" "TMe-261"  "TMe-1804" "TMe-2161" "TMe-234" 
#> [13] "TMe-3100" "TMe-3336" "TMe-1787" "TMe-1939" "TMe-3644" "TMe-123" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-1377" "TMe-1148" "TMe-3273" "TMe-3781"
#>  [7] "TMe-2971" "TMe-698"  "TMe-1434" "TMe-241"  "TMe-3428" "TMe-3255"
#> [13] "TMe-3527" "TMe-280"  "TMe-259"  "TMe-1020" "TMe-279"  "TMe-619" 
#> [19] "TMe-25"   "TMe-2567" "TMe-1988" "TMe-737"  "TMe-2956" "TMe-1297"
#> [25] "TMe-699"  "TMe-1987" "TMe-1155" "TMe-2375" "TMe-191"  "TMe-1580"
#> [31] "TMe-1144" "TMe-3072" "TMe-1376" "TMe-57"   "TMe-2318"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-723"  "TMe-419" 
#>  [7] "TMe-3329" "TMe-1004" "TMe-2355" "TMe-1880" "TMe-423"  "TMe-2853"
#> [13] "TMe-1427" "TMe-588"  "TMe-574"  "TMe-629"  "TMe-688"  "TMe-333" 
#> [19] "TMe-2124" "TMe-1375" "TMe-769"  "TMe-344"  "TMe-1220" "TMe-532" 
#> [25] "TMe-745"  "TMe-1534" "TMe-2151" "TMe-1293" "TMe-332"  "TMe-247" 
#> [31] "TMe-1012" "TMe-1290" "TMe-863"  "TMe-439"  "TMe-682"  "TMe-1131"
#> [37] "TMe-487"  "TMe-997"  "TMe-1268" "TMe-1760"
#> 
#> [[6]]
#>  [1] "TMe-1503" "TMe-1608" "TMe-2983" "TMe-1816" "TMe-693"  "TMe-1392"
#>  [7] "TMe-2791" "TMe-1646" "TMe-661"  "TMe-1124" "TMe-1566" "TMe-2543"
#> [13] "TMe-1900" "TMe-1302" "TMe-1756" "TMe-1481" "TMe-1239"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_best_sum_margalef,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-3694" "TMe-3236" "TMe-1086"
#>  [7] "TMe-865"  "TMe-2967" "TMe-3481" "TMe-3051" "TMe-2975" "TMe-1922"
#> [13] "TMe-3252" "TMe-3262" "TMe-937"  "TMe-1218" "TMe-2810" "TMe-3112"
#> [19] "TMe-2027" "TMe-778"  "TMe-1717" "TMe-3726" "TMe-1581" "TMe-2066"
#> [25] "TMe-469"  "TMe-500"  "TMe-952"  "TMe-1170" "TMe-1532" "TMe-1914"
#> [31] "TMe-1451"
#> 
#> [[2]]
#>  [1] "TMe-2329" "TMe-74"   "TMe-40"   "TMe-3530" "TMe-681"  "TMe-3766"
#>  [7] "TMe-960"  "TMe-2568" "TMe-369"  "TMe-3200" "TMe-3366" "TMe-2033"
#> [13] "TMe-3547" "TMe-2412" "TMe-2715" "TMe-3800" "TMe-2951" "TMe-196" 
#> [19] "TMe-339"  "TMe-2995" "TMe-674"  "TMe-2352" "TMe-1831" "TMe-2997"
#> [25] "TMe-796"  "TMe-1505" "TMe-1668" "TMe-1732" "TMe-3495" "TMe-2757"
#> [31] "TMe-3557"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-3133" "TMe-3397" "TMe-234"  "TMe-1804" "TMe-2161" "TMe-1738"
#> [13] "TMe-3100" "TMe-3336" "TMe-1787" "TMe-1939" "TMe-3644" "TMe-123" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-2458" "TMe-3255" "TMe-3327"
#>  [7] "TMe-1376" "TMe-698"  "TMe-1434" "TMe-280"  "TMe-3428" "TMe-3527"
#> [13] "TMe-241"  "TMe-737"  "TMe-1020" "TMe-259"  "TMe-279"  "TMe-25"  
#> [19] "TMe-619"  "TMe-3273" "TMe-1988" "TMe-1148" "TMe-2956" "TMe-1297"
#> [25] "TMe-699"  "TMe-1987" "TMe-1155" "TMe-2375" "TMe-191"  "TMe-1580"
#> [31] "TMe-1144" "TMe-3072" "TMe-57"   "TMe-2318" "TMe-1016"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-723"  "TMe-419" 
#>  [7] "TMe-3329" "TMe-1004" "TMe-2853" "TMe-1003" "TMe-423"  "TMe-1427"
#> [13] "TMe-588"  "TMe-574"  "TMe-629"  "TMe-688"  "TMe-333"  "TMe-2124"
#> [19] "TMe-1375" "TMe-769"  "TMe-344"  "TMe-1220" "TMe-532"  "TMe-745" 
#> [25] "TMe-1534" "TMe-2151" "TMe-1293" "TMe-332"  "TMe-247"  "TMe-1012"
#> [31] "TMe-1290" "TMe-863"  "TMe-439"  "TMe-682"  "TMe-1131" "TMe-487" 
#> [37] "TMe-997"  "TMe-1880" "TMe-1268" "TMe-1760"
#> 
#> [[6]]
#>  [1] "TMe-936"  "TMe-1403" "TMe-2791" "TMe-1816" "TMe-661"  "TMe-1566"
#>  [7] "TMe-2543" "TMe-1124" "TMe-3177" "TMe-693"  "TMe-751"  "TMe-1302"
#> [13] "TMe-1392" "TMe-1756" "TMe-1900" "TMe-1481" "TMe-1608"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_mean_richness,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-3694" "TMe-3236" "TMe-1086"
#>  [7] "TMe-865"  "TMe-2967" "TMe-3481" "TMe-3051" "TMe-2975" "TMe-1922"
#> [13] "TMe-3252" "TMe-3262" "TMe-937"  "TMe-1218" "TMe-2810" "TMe-3112"
#> [19] "TMe-2027" "TMe-778"  "TMe-1717" "TMe-3726" "TMe-1581" "TMe-2066"
#> [25] "TMe-469"  "TMe-500"  "TMe-952"  "TMe-1170" "TMe-1532" "TMe-1914"
#> [31] "TMe-1451"
#> 
#> [[2]]
#>  [1] "TMe-2611" "TMe-369"  "TMe-3366" "TMe-1754" "TMe-3800" "TMe-2033"
#>  [7] "TMe-3200" "TMe-40"   "TMe-674"  "TMe-796"  "TMe-74"   "TMe-2329"
#> [13] "TMe-3547" "TMe-2412" "TMe-2715" "TMe-681"  "TMe-2951" "TMe-196" 
#> [19] "TMe-339"  "TMe-2995" "TMe-2352" "TMe-1831" "TMe-2997" "TMe-2568"
#> [25] "TMe-3766" "TMe-1505" "TMe-1668" "TMe-1732" "TMe-3495" "TMe-960" 
#> [31] "TMe-2757"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-3133" "TMe-3397" "TMe-234"  "TMe-1804" "TMe-2161" "TMe-1738"
#> [13] "TMe-3100" "TMe-3336" "TMe-1787" "TMe-1939" "TMe-3644" "TMe-123" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-2807" "TMe-3255" "TMe-698" 
#>  [7] "TMe-266"  "TMe-1434" "TMe-3248" "TMe-1020" "TMe-3527" "TMe-1579"
#> [13] "TMe-241"  "TMe-737"  "TMe-279"  "TMe-25"   "TMe-3428" "TMe-619" 
#> [19] "TMe-3273" "TMe-280"  "TMe-1148" "TMe-2956" "TMe-1297" "TMe-699" 
#> [25] "TMe-259"  "TMe-1987" "TMe-1155" "TMe-2375" "TMe-191"  "TMe-1580"
#> [31] "TMe-1144" "TMe-3072" "TMe-1376" "TMe-57"   "TMe-2318"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-723"  "TMe-419" 
#>  [7] "TMe-3329" "TMe-1004" "TMe-2853" "TMe-1003" "TMe-423"  "TMe-1427"
#> [13] "TMe-588"  "TMe-574"  "TMe-629"  "TMe-688"  "TMe-333"  "TMe-2124"
#> [19] "TMe-1375" "TMe-769"  "TMe-344"  "TMe-1220" "TMe-532"  "TMe-745" 
#> [25] "TMe-1534" "TMe-2151" "TMe-1293" "TMe-332"  "TMe-247"  "TMe-1012"
#> [31] "TMe-1290" "TMe-863"  "TMe-439"  "TMe-682"  "TMe-1131" "TMe-487" 
#> [37] "TMe-997"  "TMe-1880" "TMe-1268" "TMe-1760"
#> 
#> [[6]]
#>  [1] "TMe-1992" "TMe-531"  "TMe-1661" "TMe-693"  "TMe-1392" "TMe-1816"
#>  [7] "TMe-1124" "TMe-625"  "TMe-751"  "TMe-1608" "TMe-1302" "TMe-1900"
#> [13] "TMe-1592" "TMe-2791" "TMe-1566" "TMe-1481" "TMe-1239"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_sum_richness,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2453" "TMe-2785"
#>  [7] "TMe-1096" "TMe-778"  "TMe-3112" "TMe-41"   "TMe-3130" "TMe-1086"
#> [13] "TMe-2027" "TMe-1190" "TMe-3481" "TMe-3115" "TMe-3262" "TMe-1922"
#> [19] "TMe-1218" "TMe-3252" "TMe-3437" "TMe-1869" "TMe-2237" "TMe-566" 
#> [25] "TMe-2975" "TMe-3398" "TMe-3111" "TMe-2810" "TMe-569"  "TMe-937" 
#> [31] "TMe-3236"
#> 
#> [[2]]
#>  [1] "TMe-796"  "TMe-3766" "TMe-2412" "TMe-369"  "TMe-2329" "TMe-3200"
#>  [7] "TMe-1668" "TMe-1385" "TMe-2033" "TMe-1754" "TMe-3547" "TMe-74"  
#> [13] "TMe-2952" "TMe-171"  "TMe-539"  "TMe-2568" "TMe-40"   "TMe-2352"
#> [19] "TMe-3101" "TMe-3530" "TMe-2021" "TMe-3447" "TMe-674"  "TMe-3284"
#> [25] "TMe-3805" "TMe-3366" "TMe-3093" "TMe-2257" "TMe-2211" "TMe-2000"
#> [31] "TMe-1831"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-234"  "TMe-261"  "TMe-3133" "TMe-1804" "TMe-3100" "TMe-425" 
#> [13] "TMe-3631" "TMe-2161" "TMe-1230" "TMe-1819" "TMe-2374" "TMe-3148"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-956"  "TMe-698" 
#>  [7] "TMe-1434" "TMe-3527" "TMe-2807" "TMe-241"  "TMe-3542" "TMe-1020"
#> [13] "TMe-2958" "TMe-1988" "TMe-2924" "TMe-3273" "TMe-619"  "TMe-18"  
#> [19] "TMe-2971" "TMe-3428" "TMe-812"  "TMe-1377" "TMe-280"  "TMe-266" 
#> [25] "TMe-3189" "TMe-3242" "TMe-2956" "TMe-1700" "TMe-1776" "TMe-1123"
#> [31] "TMe-25"   "TMe-1376" "TMe-279"  "TMe-427"  "TMe-3390"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-419"  "TMe-723" 
#>  [7] "TMe-3329" "TMe-2003" "TMe-2853" "TMe-877"  "TMe-2290" "TMe-423" 
#> [13] "TMe-1375" "TMe-2355" "TMe-627"  "TMe-755"  "TMe-532"  "TMe-1220"
#> [19] "TMe-2907" "TMe-1541" "TMe-167"  "TMe-1004" "TMe-2425" "TMe-1500"
#> [25] "TMe-816"  "TMe-1427" "TMe-2750" "TMe-730"  "TMe-667"  "TMe-487" 
#> [31] "TMe-1534" "TMe-344"  "TMe-98"   "TMe-1196" "TMe-1934" "TMe-1290"
#> [37] "TMe-823"  "TMe-2855" "TMe-2192" "TMe-439" 
#> 
#> [[6]]
#>  [1] "TMe-310"  "TMe-2791" "TMe-1992" "TMe-1756" "TMe-2543" "TMe-751" 
#>  [7] "TMe-1592" "TMe-693"  "TMe-1816" "TMe-1566" "TMe-1124" "TMe-1392"
#> [13] "TMe-1403" "TMe-2818" "TMe-2983" "TMe-1661" "TMe-1035"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_mean_shannon,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2453" "TMe-2785"
#>  [7] "TMe-1096" "TMe-778"  "TMe-3112" "TMe-41"   "TMe-3130" "TMe-1086"
#> [13] "TMe-2027" "TMe-1190" "TMe-3481" "TMe-3115" "TMe-3262" "TMe-1922"
#> [19] "TMe-1218" "TMe-3252" "TMe-3437" "TMe-1869" "TMe-2237" "TMe-566" 
#> [25] "TMe-2975" "TMe-3398" "TMe-3111" "TMe-2810" "TMe-569"  "TMe-937" 
#> [31] "TMe-3236"
#> 
#> [[2]]
#>  [1] "TMe-3530" "TMe-1831" "TMe-3766" "TMe-40"   "TMe-2329" "TMe-1385"
#>  [7] "TMe-3200" "TMe-171"  "TMe-2952" "TMe-2033" "TMe-369"  "TMe-2611"
#> [13] "TMe-2568" "TMe-1754" "TMe-3366" "TMe-3805" "TMe-3101" "TMe-74"  
#> [19] "TMe-539"  "TMe-674"  "TMe-2021" "TMe-3547" "TMe-3447" "TMe-1668"
#> [25] "TMe-2352" "TMe-3284" "TMe-3093" "TMe-2257" "TMe-2000" "TMe-2211"
#> [31] "TMe-796" 
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-234"  "TMe-261"  "TMe-3133" "TMe-1804" "TMe-3100" "TMe-425" 
#> [13] "TMe-3631" "TMe-2161" "TMe-1230" "TMe-1819" "TMe-2374" "TMe-3148"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-956"  "TMe-698" 
#>  [7] "TMe-1434" "TMe-3527" "TMe-2807" "TMe-241"  "TMe-3542" "TMe-1020"
#> [13] "TMe-2958" "TMe-1988" "TMe-2924" "TMe-3273" "TMe-619"  "TMe-18"  
#> [19] "TMe-2971" "TMe-3428" "TMe-812"  "TMe-1377" "TMe-280"  "TMe-266" 
#> [25] "TMe-3189" "TMe-3242" "TMe-2956" "TMe-1700" "TMe-1776" "TMe-1123"
#> [31] "TMe-25"   "TMe-1376" "TMe-279"  "TMe-427"  "TMe-3390"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-419"  "TMe-723" 
#>  [7] "TMe-3329" "TMe-2003" "TMe-2853" "TMe-877"  "TMe-2290" "TMe-423" 
#> [13] "TMe-1375" "TMe-2355" "TMe-627"  "TMe-755"  "TMe-532"  "TMe-1220"
#> [19] "TMe-2907" "TMe-1541" "TMe-167"  "TMe-1004" "TMe-2425" "TMe-1500"
#> [25] "TMe-816"  "TMe-1427" "TMe-2750" "TMe-730"  "TMe-667"  "TMe-487" 
#> [31] "TMe-1534" "TMe-344"  "TMe-98"   "TMe-1196" "TMe-1934" "TMe-1290"
#> [37] "TMe-823"  "TMe-2855" "TMe-2192" "TMe-439" 
#> 
#> [[6]]
#>  [1] "TMe-1062" "TMe-1124" "TMe-2791" "TMe-1566" "TMe-1816" "TMe-1661"
#>  [7] "TMe-1302" "TMe-1392" "TMe-2983" "TMe-1403" "TMe-693"  "TMe-1035"
#> [13] "TMe-2543" "TMe-751"  "TMe-2818" "TMe-1900" "TMe-310" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_sum_shannon,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-566" 
#>  [7] "TMe-2453" "TMe-41"   "TMe-3252" "TMe-1096" "TMe-500"  "TMe-3112"
#> [13] "TMe-1190" "TMe-3437" "TMe-2237" "TMe-3115" "TMe-1922" "TMe-778" 
#> [19] "TMe-1086" "TMe-3685" "TMe-2027" "TMe-2810" "TMe-1960" "TMe-3726"
#> [25] "TMe-3262" "TMe-3398" "TMe-1869" "TMe-569"  "TMe-3236" "TMe-937" 
#> [31] "TMe-2975"
#> 
#> [[2]]
#>  [1] "TMe-3284" "TMe-2329" "TMe-1385" "TMe-2611" "TMe-1754" "TMe-369" 
#>  [7] "TMe-3766" "TMe-2033" "TMe-2952" "TMe-2568" "TMe-171"  "TMe-3200"
#> [13] "TMe-539"  "TMe-3530" "TMe-2352" "TMe-40"   "TMe-3101" "TMe-3805"
#> [19] "TMe-674"  "TMe-2021" "TMe-3447" "TMe-3366" "TMe-2903" "TMe-251" 
#> [25] "TMe-2257" "TMe-3547" "TMe-2000" "TMe-2211" "TMe-1831" "TMe-1732"
#> [31] "TMe-3258"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1863" "TMe-3148" "TMe-2161" "TMe-3631" "TMe-425"  "TMe-1804"
#> [13] "TMe-1819" "TMe-261"  "TMe-3100" "TMe-2374" "TMe-1230" "TMe-3397"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-241"  "TMe-2924" "TMe-812" 
#>  [7] "TMe-3542" "TMe-1377" "TMe-2971" "TMe-3273" "TMe-619"  "TMe-3527"
#> [13] "TMe-1434" "TMe-698"  "TMe-3242" "TMe-280"  "TMe-1020" "TMe-2958"
#> [19] "TMe-3189" "TMe-1776" "TMe-1988" "TMe-3390" "TMe-3255" "TMe-1376"
#> [25] "TMe-1297" "TMe-3198" "TMe-1350" "TMe-3428" "TMe-2956" "TMe-1144"
#> [31] "TMe-25"   "TMe-372"  "TMe-259"  "TMe-266"  "TMe-1579"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-2290"
#>  [7] "TMe-2355" "TMe-2003" "TMe-755"  "TMe-419"  "TMe-3329" "TMe-627" 
#> [13] "TMe-1375" "TMe-2853" "TMe-423"  "TMe-1541" "TMe-2750" "TMe-344" 
#> [19] "TMe-823"  "TMe-1534" "TMe-730"  "TMe-2192" "TMe-1500" "TMe-532" 
#> [25] "TMe-487"  "TMe-2907" "TMe-877"  "TMe-1196" "TMe-167"  "TMe-1159"
#> [31] "TMe-795"  "TMe-98"   "TMe-667"  "TMe-2855" "TMe-816"  "TMe-224" 
#> [37] "TMe-1934" "TMe-1265" "TMe-1427" "TMe-688" 
#> 
#> [[6]]
#>  [1] "TMe-1992" "TMe-1445" "TMe-1661" "TMe-2983" "TMe-1566" "TMe-2791"
#>  [7] "TMe-1816" "TMe-1403" "TMe-2818" "TMe-693"  "TMe-1124" "TMe-1035"
#> [13] "TMe-2543" "TMe-751"  "TMe-598"  "TMe-310"  "TMe-1392"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_mean_simpson,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-1190"
#>  [7] "TMe-2453" "TMe-41"   "TMe-1589" "TMe-3262" "TMe-1096" "TMe-3481"
#> [13] "TMe-1086" "TMe-3115" "TMe-3437" "TMe-778"  "TMe-2027" "TMe-3252"
#> [19] "TMe-566"  "TMe-3112" "TMe-1922" "TMe-1869" "TMe-2810" "TMe-3726"
#> [25] "TMe-2237" "TMe-1218" "TMe-569"  "TMe-3111" "TMe-2975" "TMe-300" 
#> [31] "TMe-3398"
#> 
#> [[2]]
#>  [1] "TMe-2021" "TMe-3766" "TMe-369"  "TMe-2568" "TMe-171"  "TMe-2952"
#>  [7] "TMe-2033" "TMe-2329" "TMe-3200" "TMe-3530" "TMe-40"   "TMe-2352"
#> [13] "TMe-1385" "TMe-3101" "TMe-1754" "TMe-539"  "TMe-674"  "TMe-3284"
#> [19] "TMe-3805" "TMe-2257" "TMe-3366" "TMe-2211" "TMe-2000" "TMe-3447"
#> [25] "TMe-3547" "TMe-251"  "TMe-2903" "TMe-1831" "TMe-3093" "TMe-1732"
#> [31] "TMe-3258"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1863" "TMe-3148" "TMe-2161" "TMe-3631" "TMe-425"  "TMe-1804"
#> [13] "TMe-2374" "TMe-3100" "TMe-261"  "TMe-1230" "TMe-1819" "TMe-3397"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-241"  "TMe-2971" "TMe-1377"
#>  [7] "TMe-2924" "TMe-3255" "TMe-698"  "TMe-1434" "TMe-3527" "TMe-1020"
#> [13] "TMe-3542" "TMe-3273" "TMe-3242" "TMe-619"  "TMe-3390" "TMe-280" 
#> [19] "TMe-812"  "TMe-3189" "TMe-2958" "TMe-1988" "TMe-1376" "TMe-1776"
#> [25] "TMe-1297" "TMe-3198" "TMe-1350" "TMe-3428" "TMe-2956" "TMe-1144"
#> [31] "TMe-25"   "TMe-372"  "TMe-259"  "TMe-266"  "TMe-1579"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-2290"
#>  [7] "TMe-2355" "TMe-2003" "TMe-755"  "TMe-419"  "TMe-3329" "TMe-627" 
#> [13] "TMe-1375" "TMe-2853" "TMe-423"  "TMe-1541" "TMe-2750" "TMe-344" 
#> [19] "TMe-823"  "TMe-1534" "TMe-730"  "TMe-2192" "TMe-1500" "TMe-532" 
#> [25] "TMe-487"  "TMe-2907" "TMe-877"  "TMe-1196" "TMe-167"  "TMe-1159"
#> [31] "TMe-795"  "TMe-98"   "TMe-667"  "TMe-2855" "TMe-816"  "TMe-224" 
#> [37] "TMe-1934" "TMe-1265" "TMe-1427" "TMe-688" 
#> 
#> [[6]]
#>  [1] "TMe-310"  "TMe-1035" "TMe-1062" "TMe-1403" "TMe-2818" "TMe-1816"
#>  [7] "TMe-751"  "TMe-2791" "TMe-1566" "TMe-2983" "TMe-2543" "TMe-693" 
#> [13] "TMe-1124" "TMe-598"  "TMe-1661" "TMe-1392" "TMe-1992"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_sum_simpson,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-2453"
#>  [7] "TMe-1190" "TMe-41"   "TMe-3685" "TMe-2027" "TMe-3719" "TMe-1096"
#> [13] "TMe-778"  "TMe-1086" "TMe-3115" "TMe-3262" "TMe-566"  "TMe-3112"
#> [19] "TMe-1922" "TMe-3437" "TMe-3539" "TMe-3236" "TMe-569"  "TMe-2237"
#> [25] "TMe-2810" "TMe-300"  "TMe-1589" "TMe-3398" "TMe-1869" "TMe-865" 
#> [31] "TMe-2975"
#> 
#> [[2]]
#>  [1] "TMe-3366" "TMe-1754" "TMe-369"  "TMe-3258" "TMe-2033" "TMe-2952"
#>  [7] "TMe-3766" "TMe-1385" "TMe-2329" "TMe-3200" "TMe-40"   "TMe-3530"
#> [13] "TMe-171"  "TMe-2568" "TMe-3805" "TMe-1668" "TMe-2211" "TMe-3101"
#> [19] "TMe-2352" "TMe-539"  "TMe-2021" "TMe-674"  "TMe-3284" "TMe-2257"
#> [25] "TMe-3547" "TMe-2000" "TMe-2903" "TMe-251"  "TMe-3447" "TMe-3093"
#> [31] "TMe-796" 
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1230" "TMe-1863" "TMe-1819" "TMe-261"  "TMe-3148" "TMe-1804"
#> [13] "TMe-425"  "TMe-3631" "TMe-2161" "TMe-2374" "TMe-3397" "TMe-3100"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-2807" "TMe-1434"
#>  [7] "TMe-3527" "TMe-698"  "TMe-2971" "TMe-241"  "TMe-3409" "TMe-1988"
#> [13] "TMe-2924" "TMe-3273" "TMe-1376" "TMe-2958" "TMe-1020" "TMe-3542"
#> [19] "TMe-3242" "TMe-280"  "TMe-812"  "TMe-1776" "TMe-3189" "TMe-2956"
#> [25] "TMe-3198" "TMe-1377" "TMe-3428" "TMe-266"  "TMe-3390" "TMe-1297"
#> [31] "TMe-1579" "TMe-619"  "TMe-25"   "TMe-1350" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-2290"
#>  [7] "TMe-2355" "TMe-2003" "TMe-755"  "TMe-419"  "TMe-3329" "TMe-627" 
#> [13] "TMe-1220" "TMe-2192" "TMe-2853" "TMe-423"  "TMe-1375" "TMe-532" 
#> [19] "TMe-877"  "TMe-1500" "TMe-2750" "TMe-730"  "TMe-167"  "TMe-2907"
#> [25] "TMe-816"  "TMe-667"  "TMe-1541" "TMe-487"  "TMe-1534" "TMe-344" 
#> [31] "TMe-1196" "TMe-920"  "TMe-98"   "TMe-1401" "TMe-1159" "TMe-823" 
#> [37] "TMe-1934" "TMe-224"  "TMe-1004" "TMe-1295"
#> 
#> [[6]]
#>  [1] "TMe-310"  "TMe-693"  "TMe-1403" "TMe-1035" "TMe-1124" "TMe-1816"
#>  [7] "TMe-2818" "TMe-1566" "TMe-2983" "TMe-1661" "TMe-2791" "TMe-751" 
#> [13] "TMe-2543" "TMe-1392" "TMe-620"  "TMe-1445" "TMe-1992"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_mean_mcintosh,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-2453"
#>  [7] "TMe-1190" "TMe-41"   "TMe-3685" "TMe-2027" "TMe-3719" "TMe-1096"
#> [13] "TMe-778"  "TMe-1086" "TMe-3115" "TMe-3262" "TMe-566"  "TMe-3112"
#> [19] "TMe-1922" "TMe-3437" "TMe-3539" "TMe-3236" "TMe-569"  "TMe-2237"
#> [25] "TMe-2810" "TMe-300"  "TMe-1589" "TMe-3398" "TMe-1869" "TMe-865" 
#> [31] "TMe-2975"
#> 
#> [[2]]
#>  [1] "TMe-2033" "TMe-369"  "TMe-3557" "TMe-2952" "TMe-1385" "TMe-1754"
#>  [7] "TMe-2611" "TMe-2329" "TMe-2568" "TMe-3766" "TMe-40"   "TMe-3530"
#> [13] "TMe-3101" "TMe-3200" "TMe-2021" "TMe-171"  "TMe-3805" "TMe-1831"
#> [19] "TMe-539"  "TMe-674"  "TMe-3366" "TMe-2257" "TMe-2352" "TMe-3284"
#> [25] "TMe-2903" "TMe-3547" "TMe-2000" "TMe-3447" "TMe-3093" "TMe-2211"
#> [31] "TMe-796" 
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1230" "TMe-1863" "TMe-1819" "TMe-261"  "TMe-3148" "TMe-1804"
#> [13] "TMe-425"  "TMe-3631" "TMe-2161" "TMe-2374" "TMe-3397" "TMe-3100"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-2807" "TMe-1434"
#>  [7] "TMe-3527" "TMe-698"  "TMe-2971" "TMe-241"  "TMe-3409" "TMe-1988"
#> [13] "TMe-2924" "TMe-3273" "TMe-1376" "TMe-2958" "TMe-1020" "TMe-3542"
#> [19] "TMe-3242" "TMe-280"  "TMe-812"  "TMe-1776" "TMe-3189" "TMe-2956"
#> [25] "TMe-3198" "TMe-1377" "TMe-3428" "TMe-266"  "TMe-3390" "TMe-1297"
#> [31] "TMe-1579" "TMe-619"  "TMe-25"   "TMe-1350" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-2290"
#>  [7] "TMe-2355" "TMe-2003" "TMe-755"  "TMe-419"  "TMe-3329" "TMe-627" 
#> [13] "TMe-1220" "TMe-2192" "TMe-2853" "TMe-423"  "TMe-1375" "TMe-532" 
#> [19] "TMe-877"  "TMe-1500" "TMe-2750" "TMe-730"  "TMe-167"  "TMe-2907"
#> [25] "TMe-816"  "TMe-667"  "TMe-1541" "TMe-487"  "TMe-1534" "TMe-344" 
#> [31] "TMe-1196" "TMe-920"  "TMe-98"   "TMe-1401" "TMe-1159" "TMe-823" 
#> [37] "TMe-1934" "TMe-224"  "TMe-1004" "TMe-1295"
#> 
#> [[6]]
#>  [1] "TMe-598"  "TMe-531"  "TMe-2983" "TMe-1816" "TMe-693"  "TMe-1403"
#>  [7] "TMe-1392" "TMe-1124" "TMe-2791" "TMe-1566" "TMe-2543" "TMe-751" 
#> [13] "TMe-1035" "TMe-310"  "TMe-2818" "TMe-1661" "TMe-1992"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_sum_mcintosh,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-2453"
#>  [7] "TMe-1190" "TMe-41"   "TMe-3481" "TMe-2027" "TMe-3130" "TMe-1096"
#> [13] "TMe-778"  "TMe-1086" "TMe-3115" "TMe-3112" "TMe-1922" "TMe-3437"
#> [19] "TMe-566"  "TMe-3262" "TMe-2237" "TMe-3539" "TMe-2810" "TMe-1869"
#> [25] "TMe-3398" "TMe-3252" "TMe-2975" "TMe-3111" "TMe-569"  "TMe-1218"
#> [31] "TMe-300" 
#> 
#> [[2]]
#>  [1] "TMe-539"  "TMe-3766" "TMe-2903" "TMe-2329" "TMe-251"  "TMe-2568"
#>  [7] "TMe-171"  "TMe-369"  "TMe-2033" "TMe-2952" "TMe-2352" "TMe-40"  
#> [13] "TMe-1385" "TMe-3200" "TMe-1754" "TMe-3530" "TMe-3101" "TMe-3547"
#> [19] "TMe-1668" "TMe-674"  "TMe-3805" "TMe-2257" "TMe-3366" "TMe-3447"
#> [25] "TMe-74"   "TMe-2021" "TMe-3093" "TMe-2000" "TMe-2211" "TMe-3258"
#> [31] "TMe-796" 
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1863" "TMe-261"  "TMe-3133" "TMe-3100" "TMe-1804" "TMe-425" 
#> [13] "TMe-3631" "TMe-1819" "TMe-1230" "TMe-2161" "TMe-2374" "TMe-3148"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-1434" "TMe-1700"
#>  [7] "TMe-698"  "TMe-3527" "TMe-2971" "TMe-241"  "TMe-1155" "TMe-2924"
#> [13] "TMe-1988" "TMe-3273" "TMe-1020" "TMe-2958" "TMe-18"   "TMe-372" 
#> [19] "TMe-280"  "TMe-812"  "TMe-3542" "TMe-1776" "TMe-3428" "TMe-1376"
#> [25] "TMe-956"  "TMe-3189" "TMe-3242" "TMe-1297" "TMe-2956" "TMe-266" 
#> [31] "TMe-25"   "TMe-1123" "TMe-3390" "TMe-1579" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-3329"
#>  [7] "TMe-419"  "TMe-2853" "TMe-2003" "TMe-755"  "TMe-423"  "TMe-2290"
#> [13] "TMe-1375" "TMe-2355" "TMe-627"  "TMe-1220" "TMe-532"  "TMe-1500"
#> [19] "TMe-877"  "TMe-2750" "TMe-2907" "TMe-167"  "TMe-816"  "TMe-1004"
#> [25] "TMe-1427" "TMe-487"  "TMe-667"  "TMe-730"  "TMe-1534" "TMe-1196"
#> [31] "TMe-344"  "TMe-98"   "TMe-2855" "TMe-823"  "TMe-688"  "TMe-2192"
#> [37] "TMe-362"  "TMe-1290" "TMe-1541" "TMe-1934"
#> 
#> [[6]]
#>  [1] "TMe-2791" "TMe-1992" "TMe-1403" "TMe-1566" "TMe-1816" "TMe-2818"
#>  [7] "TMe-751"  "TMe-693"  "TMe-2543" "TMe-1035" "TMe-1124" "TMe-1661"
#> [13] "TMe-1392" "TMe-2983" "TMe-310"  "TMe-598"  "TMe-1445"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_mean_brillouin,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-2967" "TMe-2785" "TMe-2453"
#>  [7] "TMe-1190" "TMe-41"   "TMe-3481" "TMe-2027" "TMe-3130" "TMe-1096"
#> [13] "TMe-778"  "TMe-1086" "TMe-3115" "TMe-3112" "TMe-1922" "TMe-3437"
#> [19] "TMe-566"  "TMe-3262" "TMe-2237" "TMe-3539" "TMe-2810" "TMe-1869"
#> [25] "TMe-3398" "TMe-3252" "TMe-2975" "TMe-3111" "TMe-569"  "TMe-1218"
#> [31] "TMe-300" 
#> 
#> [[2]]
#>  [1] "TMe-2257" "TMe-1385" "TMe-2952" "TMe-3766" "TMe-369"  "TMe-2033"
#>  [7] "TMe-1754" "TMe-2568" "TMe-3284" "TMe-2329" "TMe-3200" "TMe-40"  
#> [13] "TMe-3530" "TMe-3366" "TMe-3805" "TMe-171"  "TMe-1831" "TMe-3547"
#> [19] "TMe-2211" "TMe-3101" "TMe-2352" "TMe-2903" "TMe-539"  "TMe-674" 
#> [25] "TMe-2021" "TMe-74"   "TMe-3447" "TMe-3093" "TMe-2000" "TMe-3258"
#> [31] "TMe-796" 
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-234" 
#>  [7] "TMe-1863" "TMe-261"  "TMe-3133" "TMe-3100" "TMe-1804" "TMe-425" 
#> [13] "TMe-3631" "TMe-1819" "TMe-1230" "TMe-2161" "TMe-2374" "TMe-3148"
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-3255" "TMe-1434" "TMe-1700"
#>  [7] "TMe-698"  "TMe-3527" "TMe-2971" "TMe-241"  "TMe-1155" "TMe-2924"
#> [13] "TMe-1988" "TMe-3273" "TMe-1020" "TMe-2958" "TMe-18"   "TMe-372" 
#> [19] "TMe-280"  "TMe-812"  "TMe-3542" "TMe-1776" "TMe-3428" "TMe-1376"
#> [25] "TMe-956"  "TMe-3189" "TMe-3242" "TMe-1297" "TMe-2956" "TMe-266" 
#> [31] "TMe-25"   "TMe-1123" "TMe-3390" "TMe-1579" "TMe-897" 
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-2425" "TMe-723"  "TMe-3329"
#>  [7] "TMe-419"  "TMe-2853" "TMe-2003" "TMe-755"  "TMe-423"  "TMe-2290"
#> [13] "TMe-1375" "TMe-2355" "TMe-627"  "TMe-1220" "TMe-532"  "TMe-1500"
#> [19] "TMe-877"  "TMe-2750" "TMe-2907" "TMe-167"  "TMe-816"  "TMe-1004"
#> [25] "TMe-1427" "TMe-487"  "TMe-667"  "TMe-730"  "TMe-1534" "TMe-1196"
#> [31] "TMe-344"  "TMe-98"   "TMe-2855" "TMe-823"  "TMe-688"  "TMe-2192"
#> [37] "TMe-362"  "TMe-1290" "TMe-1541" "TMe-1934"
#> 
#> [[6]]
#>  [1] "TMe-1900" "TMe-1302" "TMe-1403" "TMe-1816" "TMe-2791" "TMe-1062"
#>  [7] "TMe-693"  "TMe-2543" "TMe-751"  "TMe-1661" "TMe-1035" "TMe-1124"
#> [13] "TMe-1392" "TMe-1566" "TMe-2983" "TMe-2818" "TMe-310" 
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_sum_brillouin,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-3694" "TMe-3236" "TMe-1086"
#>  [7] "TMe-865"  "TMe-41"   "TMe-3051" "TMe-2785" "TMe-2975" "TMe-1922"
#> [13] "TMe-3252" "TMe-3262" "TMe-937"  "TMe-1218" "TMe-3112" "TMe-2027"
#> [19] "TMe-778"  "TMe-1717" "TMe-3726" "TMe-1581" "TMe-2066" "TMe-469" 
#> [25] "TMe-500"  "TMe-952"  "TMe-1170" "TMe-1532" "TMe-1914" "TMe-1451"
#> [31] "TMe-1117"
#> 
#> [[2]]
#>  [1] "TMe-2033" "TMe-369"  "TMe-196"  "TMe-2952" "TMe-2568" "TMe-2021"
#>  [7] "TMe-40"   "TMe-3800" "TMe-796"  "TMe-74"   "TMe-2715" "TMe-3366"
#> [13] "TMe-3547" "TMe-3200" "TMe-2412" "TMe-681"  "TMe-2951" "TMe-339" 
#> [19] "TMe-2995" "TMe-674"  "TMe-2352" "TMe-1831" "TMe-2997" "TMe-3766"
#> [25] "TMe-1505" "TMe-2329" "TMe-1668" "TMe-1732" "TMe-3495" "TMe-960" 
#> [31] "TMe-2757"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-3133" "TMe-3397" "TMe-234"  "TMe-1804" "TMe-2161" "TMe-1738"
#> [13] "TMe-3336" "TMe-3100" "TMe-1787" "TMe-1939" "TMe-3644" "TMe-123" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-2567" "TMe-2458" "TMe-3255" "TMe-3327"
#>  [7] "TMe-1376" "TMe-698"  "TMe-1434" "TMe-280"  "TMe-3428" "TMe-3527"
#> [13] "TMe-241"  "TMe-737"  "TMe-1020" "TMe-3273" "TMe-259"  "TMe-279" 
#> [19] "TMe-25"   "TMe-619"  "TMe-1988" "TMe-1148" "TMe-2956" "TMe-1297"
#> [25] "TMe-699"  "TMe-1987" "TMe-1155" "TMe-2375" "TMe-191"  "TMe-1580"
#> [31] "TMe-1144" "TMe-3072" "TMe-57"   "TMe-2318" "TMe-1016"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-723"  "TMe-419" 
#>  [7] "TMe-3329" "TMe-1004" "TMe-2853" "TMe-1003" "TMe-423"  "TMe-1427"
#> [13] "TMe-588"  "TMe-574"  "TMe-629"  "TMe-688"  "TMe-333"  "TMe-2124"
#> [19] "TMe-1375" "TMe-769"  "TMe-344"  "TMe-1220" "TMe-532"  "TMe-745" 
#> [25] "TMe-1534" "TMe-2151" "TMe-1293" "TMe-332"  "TMe-247"  "TMe-1012"
#> [31] "TMe-1290" "TMe-863"  "TMe-439"  "TMe-682"  "TMe-1131" "TMe-487" 
#> [37] "TMe-997"  "TMe-1880" "TMe-1268" "TMe-1760"
#> 
#> [[6]]
#>  [1] "TMe-631"  "TMe-1062" "TMe-2791" "TMe-1744" "TMe-2543" "TMe-1816"
#>  [7] "TMe-751"  "TMe-1392" "TMe-693"  "TMe-1756" "TMe-1900" "TMe-1124"
#> [13] "TMe-1302" "TMe-1592" "TMe-1566" "TMe-1481" "TMe-1608"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_mean_margalef,
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
#>  [1] "TMe-1830" "TMe-2934" "TMe-1425" "TMe-3694" "TMe-3236" "TMe-1086"
#>  [7] "TMe-3481" "TMe-865"  "TMe-2967" "TMe-3051" "TMe-2975" "TMe-1922"
#> [13] "TMe-3252" "TMe-3262" "TMe-2810" "TMe-937"  "TMe-1218" "TMe-2027"
#> [19] "TMe-778"  "TMe-3112" "TMe-1717" "TMe-3726" "TMe-1581" "TMe-2066"
#> [25] "TMe-469"  "TMe-500"  "TMe-952"  "TMe-1170" "TMe-1532" "TMe-1914"
#> [31] "TMe-1451"
#> 
#> [[2]]
#>  [1] "TMe-369"  "TMe-2352" "TMe-339"  "TMe-3366" "TMe-3093" "TMe-1474"
#>  [7] "TMe-3200" "TMe-74"   "TMe-674"  "TMe-2033" "TMe-796"  "TMe-40"  
#> [13] "TMe-2329" "TMe-3547" "TMe-2412" "TMe-2715" "TMe-3800" "TMe-681" 
#> [19] "TMe-2951" "TMe-196"  "TMe-2995" "TMe-1831" "TMe-2997" "TMe-2568"
#> [25] "TMe-3766" "TMe-1505" "TMe-1668" "TMe-1732" "TMe-3495" "TMe-960" 
#> [31] "TMe-2757"
#> 
#> [[3]]
#>  [1] "TMe-1790" "TMe-3750" "TMe-3715" "TMe-3556" "TMe-70"   "TMe-1863"
#>  [7] "TMe-3133" "TMe-3397" "TMe-261"  "TMe-1804" "TMe-2161" "TMe-234" 
#> [13] "TMe-3100" "TMe-3336" "TMe-1787" "TMe-1939" "TMe-3644" "TMe-123" 
#> 
#> [[4]]
#>  [1] "TMe-801"  "TMe-3191" "TMe-1377" "TMe-1148" "TMe-3273" "TMe-3781"
#>  [7] "TMe-2971" "TMe-698"  "TMe-1434" "TMe-241"  "TMe-3428" "TMe-3255"
#> [13] "TMe-3527" "TMe-280"  "TMe-259"  "TMe-1020" "TMe-279"  "TMe-619" 
#> [19] "TMe-25"   "TMe-2567" "TMe-1988" "TMe-737"  "TMe-2956" "TMe-1297"
#> [25] "TMe-699"  "TMe-1987" "TMe-1155" "TMe-2375" "TMe-191"  "TMe-1580"
#> [31] "TMe-1144" "TMe-3072" "TMe-1376" "TMe-57"   "TMe-2318"
#> 
#> [[5]]
#>  [1] "TMe-2018" "TMe-712"  "TMe-256"  "TMe-707"  "TMe-723"  "TMe-419" 
#>  [7] "TMe-3329" "TMe-1004" "TMe-2355" "TMe-1880" "TMe-423"  "TMe-2853"
#> [13] "TMe-1427" "TMe-588"  "TMe-574"  "TMe-629"  "TMe-688"  "TMe-333" 
#> [19] "TMe-2124" "TMe-1375" "TMe-769"  "TMe-344"  "TMe-1220" "TMe-532" 
#> [25] "TMe-745"  "TMe-1534" "TMe-2151" "TMe-1293" "TMe-332"  "TMe-247" 
#> [31] "TMe-1012" "TMe-1290" "TMe-863"  "TMe-439"  "TMe-682"  "TMe-1131"
#> [37] "TMe-487"  "TMe-997"  "TMe-1268" "TMe-1760"
#> 
#> [[6]]
#>  [1] "TMe-1035" "TMe-310"  "TMe-1239" "TMe-2543" "TMe-2791" "TMe-620" 
#>  [7] "TMe-1124" "TMe-693"  "TMe-1392" "TMe-1756" "TMe-1900" "TMe-1566"
#> [13] "TMe-1592" "TMe-1302" "TMe-751"  "TMe-1481" "TMe-1608"
#> 

plot_dist(d = dist_matrix, method = "isomds",
          gp = gp_vec,
          highlight =  unlist(greedysel_first_sum_margalef,
                              use.names = FALSE)) +
  labs(title = "Greed search | 1-opt first improvement",
       subtitle = "Pooled Margalef's richness index")

```
