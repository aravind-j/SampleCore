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
  quantitative,
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

- quantitative:

  Name of columns with the quantitative traits as a character vector.

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
there are no accessions specified in `always.selected` present in the
particular cluster/group ) (Nemhauser et al. 1978; Fisher et al. 1978;
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
