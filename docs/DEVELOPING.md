# Developing

### Installation Requirements
We develop on Mac, and install the following dependencies:

```bash
# Mac dependencies for building PDFs and documentation
brew install qpdf
brew cask install basictex
brew install pandoc
# for building PDFs, ensure TexLive Manager contains the required packages
sudo tlmgr update --self
sudo tlmgr update --all
sudo tlmgr install titling framed inconsolata
sudo tlmgr install collection-fontsrecommended
```

Several R packages are also required for building:

```R
R> install.packages(c('devtools', 'knitr', 'rmarkdown', 'roxygen2', 'testthat'))
```

### Building `betabinom`:
For fully building and testing betabinom, we recommend the following:

```R
devtools::build_vignettes()
devtools::install()
devtools::build(manual=TRUE)
devtools::check(manual=TRUE,cran=TRUE)
```

### Contributing
After modifying any code, please run `Rscript scripts/ci.R` to run formatting and linting.
