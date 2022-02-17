# networker: Reproducible PPI Network Creation and Visualization in R

A suite of functions to make it simple to construct PPI networks inside of R,
with an emphasis on usability and reproducibility.

## Installation
You can install `networker` with the following (assuming you have `devtools`
already installed):
```r
devtools::install_github("https://github.com/travis-m-blimkie/networker")
```

If you run into installation issues relating to the package `biomaRt` (which is
a dependency of `networker`), you can install it with the following:
```r
if (!require("BiocManager")) { # Optionally install BiocManager if needed
  install.packages("BiocManager")
}
BiocManager::install("biomaRt")
```

## Versioning
This package makes use of [SemVer](https://semver.org/).

## Authors
Travis Blimkie is the originator and principal contributor. You can check the
list of all contributors [here](https://github.com/travis-m-blimkie/networker/graphs/contributors).

## License
This project is written under the GPLv3 license, available
[here.](https://github.com/travis-m-blimkie/networker/blob/main/LICENSE.md)
