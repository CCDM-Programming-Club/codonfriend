# codonfriend

[![Travis-CI Build Status](https://travis-ci.org/CCDM-Programming-Club/codonfriend.svg?branch=master)](https://travis-ci.org/CCDM-Programming-Club/codonfriend)
[![Coverage Status](https://img.shields.io/codecov/c/github/CCDM-Programming-Club/codonfriend/master.svg)](https://codecov.io/github/CCDM-Programming-Club/codonfriend?branch=master)

Deviations from normal codon usage may indicate some compartmentalisation of mutational load (e.g. RIP) or presence of foreign DNA (e.g. horizontal transfer). Codonfriend is a small library for working with codon data and facilitating analyses of codon usage.

## Quick Install

This is a working project so for now I suggest you clone the repository and use [`devtools`](https://github.com/r-lib/devtools) or rstudio to work with it.
You should also be able to install using `devtools`.

```r
devtools::install_github("https://github.com/CCDM-Programming-Club/codonfriend.git")
```
## Roadmap

- [ ] A way to calculate CUT and CAIs for input fasta files.
- [x] Parse existing tables from EMBOSS
- [ ] Write our own version
- [ ] Plot CAI/CUT distributions (Histograms etc.)
- [ ] Basic statistical properties of CAI distributions (e.g. quartiles, fitting a model etc)
- [ ] Define what an “interesting” gene is wrt CAI or CUT e.g 5th percentile
- [ ] Validation - where are the known effectors and LGT candidates?
- [ ] Compare codon usage between genomes.
- [ ] Clustering 3-mers of “interesting” genes, are there differences between LGT candidates and RIPd candidates.
- [ ] Identify RIP-like codon changes.
- [ ] Regional biases in codon usage or CAI.

## Contributing


Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md).
By participating in this project you agree to abide by its terms.