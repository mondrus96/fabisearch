#R package for FaBiSearch

R implementation of "FaBiSearch" for change point detection, also available on CRAN at https://CRAN.R-project.org/package=fabisearch

To install and load in R:
```
install.packages("fabisearch")
library(fabisearch)
```

A trivial example of the main change point detection function, `detect.cps()`:

```
set.seed(123)
detect.cps(sim2, rank = 3, mindist = 99, nruns = 2, nreps = 2)
```

To cite package ‘fabisearch’ in publications use:

  Ondrus M, Cribben I (2023). _fabisearch: Change Point Detection in High-Dimensional Time Series Networks_. R package version 0.0.4.5,
  <https://CRAN.R-project.org/package=fabisearch>.

A BibTeX entry for LaTeX users is

```
@Manual{,
  title = {fabisearch: Change Point Detection in High-Dimensional Time Series Networks},
  author = {Martin Ondrus and Ivor Cribben},
  year = {2023},
  note = {R package version 0.0.4.5},
  url = {https://CRAN.R-project.org/package=fabisearch},
}
```
