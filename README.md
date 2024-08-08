# betabinomial

## About
This is a maintained fork of [countdata](https://cran.r-project.org/web/packages/countdata/index.html). 
It was originally published by [Thang Pham](https://orcid.org/0000-0003-0333-2492) (t.pham@amsterdamumc.nl) under a [BSD 3-Clause](./metadata/BSD_3_clause) [license](./metadata/LICENSE).

## Docs
- [Reference manual](https://cran.r-project.org/web/packages/countdata/countdata.pdf)
- [Vignettes](https://cran.r-project.org/web/packages/countdata/vignettes/countdata.html)

## Dependencies
- rmarkdown
- knitr

## Development
Install requirements:
```console
brew install qpdf
brew cask install basictex
brew install pandoc
sudo tlmgr update --self
sudo tlmgr update --all
sudo tlmgr install titling framed inconsolata
sudo tlmgr install collection-fontsrecommended
```

Building `betabinom`:
```R
install.packages('devtools')
devtools::build_vignettes()
devtools::install('betabinom')
library('betabinom')
```

