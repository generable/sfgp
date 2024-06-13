# The 'sfgp' R package 

SF + GP modeling.

## Authentication

You need to set the `GITHUB_PAT` environment variable, which should be your GitHub [personal
access token](https://github.com/settings/tokens?type=beta). 

##  Installation

* Install `cmdstanr` following the instructions [here](https://mc-stan.org/cmdstanr/).
* Install `sfgp` using

```r
remotes::install_github("generable/sfgp", ref = "main", build_vignettes = TRUE)
```

## Getting started

 More info is in documentation, that you can view with
```r
library(sfgp)
?sfgp
```
