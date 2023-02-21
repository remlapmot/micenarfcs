# A fork of a fork

This is a fork of the <https://github.com/moreno-betancur/mice> repository which enables NARFCS in the **mice** package.
That repository is itself a fork of the [**mice**](https://github.com/amices/mice/) repository.

## Why do this?

This fork renames the package to **micenarfcs**. This means that a user can have both the official **mice** package and this package installed within R at the same time (and that this package no longer overwrites the official version of **mice** on installation).

## Installation

``` r
# install.packages("remotes") # uncomment on first run
remotes::install_github("remlapmot/micenarfcs")
```
