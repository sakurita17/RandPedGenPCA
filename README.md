

# Randomized Pedigree Principal Component Analysis

A pedigree can be represented as a genetic relationship matrix (GRM). Therefore,
pedigrees can in principle be subjected to (and visualised by) principal
component analysis (PCA). However, doing this naively is very slow for large
pedigrees of, say, one million individuals. We present the Randomized Pedigree Principal
Component Analysis approach, which performs rapid PCA of a pedigree
GRM using randomized linear algebra.

[Henderson (1975)](https://doi.org/10.3168/jds.S0022-0302\(75\)84776-X) developed an
efficient way to compute the lower Cholesky factor of the inverse GRM and
[Colleau (2002)](https://doi.org/10.1186/1297-9686-34-4-409) explicitly showed
how to multiply pedigree GRM with an arbitrary vector efficiently. Our approach
uses these two algorithmic ingredients to rapidly compute the principal
components of pedigree GRM without forming the pedigree GRM. This is achieved
via the randomized singular value decomposition (rSVD) described in [Halko et
al. (2011)](http://dx.doi.org/10.1137/090771806). The resulting principal
components can reveal the underlying population structure of a pedigree.

# Paper
**Lee H, Craddock RF, Gorjanc G & Becher H (2025)**  
randPedPCA: Rapid approximation of principal components from large pedigrees  
_Genetics Selection Evolution_, 57(1), 46 [here](https://doi.org/10.1186/s12711-025-00994-y)

# R package

## Setup
`randPedPCA` is on CRAN!

[![](https://cranlogs.r-pkg.org/badges/randPedPCA)](https://cran.r-project.org/package=randPedPCA)
```
install.package("randPedPCA")
```
## First steps
For a demonstration, check out
```
library(randPedPCA)
vignette("pedigree-pca")
```

# Python prototype

The initial prototype was developed in Python. An example can be found in `notebook/Example.ipynb`, which uses the `rppca` module in this repository.

# Repo structure
```
datasets   ... datasets for testing and code for generating data (irrelavant to users)
notebooks  ... demonstration of the Python prototype
randPedPCA ... the R PACKAGE
rppca      ... Python prototype
```
