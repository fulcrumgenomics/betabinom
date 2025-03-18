# betabinom: The BetaBinomial Test
<!-- badges: start -->
[![R-CMD-check](https://github.com/fulcrumgenomics/betabinomial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fulcrumgenomics/betabinomial/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/fulcrumgenomics/betabinomial/graph/badge.svg)](https://app.codecov.io/gh/fulcrumgenomics/betabinomial)
<!-- badges: end -->

<p>
<a href float="left"="https://fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="100"/></a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with betabinom and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>


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
