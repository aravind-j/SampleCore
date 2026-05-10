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
allocation provided. Entries listed as `always.selected` are mandatorily
included in the selection. Warnings are issued if requested allocation
is smaller than the number of always-selected entries in a cluster/group
and/or when the cluster/group does not contain enough remaining entries
to fulfill the allocation.
