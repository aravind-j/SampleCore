# Version 0.1.0 - Second submission

> `DESCRIPTION`: Adding references to description field of DESCRIPTION

As there are several diverse methods implemented, I would prefer to keep the references to the documentation and mention that in the description field of DESCRIPTION.

> `R/select.distance.R`: Please make sure that you do not change the user's options, par or working directory.

Removed the par settings/option change line as it was stale code as I had moved the latter plotting with `plot_dist` to `ggplot2` with facets. 

> `R/select.diversity.R`: You write information messages to the console that cannot be easily
suppressed.

Replaced print() with message().

### Test environments
* local Windows 10 Pro 25H2, R-release (R 4.6.0) & R-devel (R 4.7.0 Pre-release).
* local Ubuntu 20.04, R-release (R 4.6.0) & R-devel (R 4.7.0 Pre-release).
* win-builder, R-release (R 4.6.0) & R-devel (R 4.7.0 Pre-release).
* github macOS Sequoia 15.7.4, R-release (R 4.6.0).
* github Ubuntu 24.04.4, R-release (R 4.6.0), R-oldrel-1 (R 4.5.3) & R-devel (R 4.7.0 Pre-release).

### R CMD check results
* There were no NOTES, ERRORs or WARNINGs.
