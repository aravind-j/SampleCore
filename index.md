## `SampleCore`: Sampling Strategies for Constructing Core Collections

![logo](https://raw.githubusercontent.com/aravind-j/SampleCore/master/inst/extdata/SampleCore.png)

###### Version : [0.1.0.9000](https://aravind-j.github.io/SampleCore/); License: [GPL-2\|GPL-3](https://www.r-project.org/Licenses/)

##### *Aravind, J., Roy, Suman and Singh, Anju M.*

Division of Germplasm Conservation, ICAR-National Bureau of Plant
Genetic Resources, New Delhi.

------------------------------------------------------------------------

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg?logo=R)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/SampleCore)](https://cran.r-project.org/package=SampleCore)
[![Dependencies](https://tinyverse.netlify.app/status/SampleCore)](https://cran.r-project.org/package=SampleCore)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/SampleCore?color=green)](https://CRAN.R-project.org/package=SampleCore)
[![develVersion](https://img.shields.io/badge/devel%20version-0.1.0.9000-orange.svg)](https://github.com/aravind-j/SampleCore)
[![Github Code
Size](https://img.shields.io/github/languages/code-size/aravind-j/SampleCore.svg)](https://github.com/aravind-j/SampleCore)
[![R-CMD-check](https://github.com/aravind-j/SampleCore/workflows/R-CMD-check/badge.svg)](https://github.com/aravind-j/SampleCore/actions)
[![Project Status:
WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![lifecycle](https://lifecycle.r-lib.org/articles/figures/lifecycle-maturing.svg)](https://lifecycle.r-lib.org/articles/stages.html#maturing)
[![Last-changedate](https://img.shields.io/badge/last%20change-2026--05--29-yellowgreen.svg)](https://github.com/aravind-j/SampleCore/)
[![Zenodo
DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20449499.svg)](https://doi.org/10.5281/zenodo.20449499)
[![Website -
pkgdown](https://img.shields.io/website-up-down-green-red/https/aravind-j.github.io/SampleCore.svg)](https://aravind-j.github.io/SampleCore/)
[![GoatCounter](https://samplecore-gh.goatcounter.com/count?p=/test)](https://samplecore-gh.goatcounter.com/)

------------------------------------------------------------------------

## Description

Implements multiple allocation and selection strategies of sampling to
construct core collections primarily from clustered or grouped germplasm
collection data. Provides methods for allocating entries to
clusters/groups based on group sizes, group-wise distance-based
diversity metrics, and group-wise diversity index estimates. Includes
procedures for selecting entries within clusters/groups through random
sampling, genetic distance-based approaches, and optimized diversity
metric–based selection methods. See the package documentation for more,
including full list of references for the methods implemented.

## Installation

The package can be installed from CRAN as follows:

The development version can be installed from github as follows:

``` r

# Install development version from Github
devtools::install_github("aravind-j/SampleCore")
```

## What’s new

To know whats new in this version type:

``` r

news(package='SampleCore')
```

## Links

[CRAN page](https://cran.r-project.org/package=SampleCore)

[Github page](https://github.com/aravind-j/SampleCore)

[Documentation website](https://aravind-j.github.io/SampleCore/)

[Zenodo DOI](https://doi.org/10.5281/zenodo.20449499)

## CRAN checks

[![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)](https://cran.r-project.org/web/checks/check_results_SampleCore.html)

| Flavour | CRAN check |
|----|----|
| r-devel-linux-x86_64-debian-clang | [![CRAN check - r-devel-linux-x86_64-debian-clang](https://badges.cranchecks.info/flavor/r-devel-linux-x86_64-debian-clang/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |
| r-devel-linux-x86_64-debian-gcc | [![CRAN check - r-devel-linux-x86_64-debian-gcc](https://badges.cranchecks.info/flavor/r-devel-linux-x86_64-debian-gcc/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |
| r-devel-linux-x86_64-fedora-clang | [![CRAN check - r-devel-linux-x86_64-fedora-clang](https://badges.cranchecks.info/flavor/r-devel-linux-x86_64-fedora-clang/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |
| r-devel-linux-x86_64-fedora-gcc | [![CRAN check - r-devel-linux-x86_64-fedora-gcc](https://badges.cranchecks.info/flavor/r-devel-linux-x86_64-fedora-gcc/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |
| r-patched-linux-x86_64 | [![CRAN check - r-patched-linux-x86_64](https://badges.cranchecks.info/flavor/r-patched-linux-x86_64/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |
| r-release-linux-x86_64 | [![CRAN check - r-release-linux-x86_64](https://badges.cranchecks.info/flavor/r-release-linux-x86_64/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |

[![Windows](https://img.shields.io/badge/Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white)](https://cran.r-project.org/web/checks/check_results_SampleCore.html)

| Flavour | CRAN check |
|----|----|
| r-devel-windows-x86_64 | [![CRAN check - r-devel-windows-x86_64](https://badges.cranchecks.info/flavor/r-devel-windows-x86_64/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |
| r-release-windows-x86_64 | [![CRAN check - r-release-windows-x86_64](https://badges.cranchecks.info/flavor/r-release-windows-x86_64/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |
| r-oldrel-windows-x86_64 | [![CRAN check - r-oldrel-windows-x86_64](https://badges.cranchecks.info/flavor/r-oldrel-windows-x86_64/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |

[![MacOS](https://img.shields.io/badge/mac%20os-000000?style=for-the-badge&logo=apple&logoColor=white)](https://cran.r-project.org/web/checks/check_results_SampleCore.html)

| Flavour | CRAN check |
|----|----|
| r-release-macos-x86_64 | [![CRAN check - r-release-macos-x86_64](https://badges.cranchecks.info/flavor/r-release-macos-x86_64/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |
| r-oldrel-macos-x86_64 | [![CRAN check - r-oldrel-macos-x86_64](https://badges.cranchecks.info/flavor/r-oldrel-macos-x86_64/SampleCore.svg)](https://cran.r-project.org/web/checks/check_results_SampleCore.html) |

## Citing `SampleCore`

To cite the methods in the package use:

``` r

citation("SampleCore")
```

``` R
To cite the R package 'SampleCore' in publications use:

  Aravind, J., Roy, S., and Singh, A. M. ().  SampleCore: Sampling Strategies for Constructing Core
  Collections. R package version 0.1.0.9000, https://aravind-j.github.io/SampleCore/.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {SampleCore: Sampling Strategies for Constructing Core Collections},
    author = {J. Aravind and Suman Roy and Anju Mahendru Singh},
    note = {R package version 0.1.0.9000 https://aravind-j.github.io/SampleCore/},
  }

This free and open-source software implements academic research by the authors and co-workers. If you use
it, please support the project by citing the package.
```
