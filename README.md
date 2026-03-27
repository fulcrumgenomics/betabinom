# betabinom: The BetaBinomial Test
<!-- badges: start -->
[![R-CMD-check](https://github.com/fulcrumgenomics/betabinomial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fulcrumgenomics/betabinomial/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/fulcrumgenomics/betabinomial/graph/badge.svg)](https://app.codecov.io/gh/fulcrumgenomics/betabinomial)
<!-- badges: end -->

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset=".github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset=".github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src=".github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with betabinom and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

## About
This is a maintained fork of [countdata](https://CRAN.R-project.org/package=countdata). 
It was originally published by [Thang Pham](https://orcid.org/0000-0003-0333-2492) (t.pham@amsterdamumc.nl) under a [BSD 3-Clause](https://opensource.org/license/bsd-3-clause) [license](LICENSE).

## Installation
### 1. Github Source
First, download the source from Github and start an `R` session:

```bash
git clone https://github.com/fulcrumgenomics/betabinom
cd betabinom
R
```

Then, download the dependencies and install betabinom:

```R
install.packages(c('devtools', 'knitr', 'rmarkdown', 'roxygen2', 'testthat'))
devtools::install()
```

### 2. CRAN
We eventually plan to add this project to CRAN.

### 3. conda-forge
We eventually plan to add this project to conda-forge.

## Docs
- [Reference manual](https://CRAN.R-project.org/package=countdata/countdata.pdf)
- [Vignettes](https://CRAN.R-project.org/package=countdata/vignettes/countdata.html)

## Development
For information on developing, please see [here](docs/DEVELOPING.md).
