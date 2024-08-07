# FOLD
This is the GitHub repository for Bayesian Clustering via Fusing of Localized Densities. In this article, we present Fusing of Localized Densities (FOLD), a decision theoretic clustering method that focuses on grouping the _localized densities_ from a mixture model. Here, you will find the ```foldcluster``` R package and code required to reproduce all plots, the simulation studies, and the application to the cell line dataset.

## R Package
To install our R package, run the following code in your console.
```r
install.packages("devtools")
devtools::install_github("adombowsky/FOLD/foldcluster", build_vignettes=TRUE)
```

## Introductory Vignette
We recommend first running the vignette ```introduction```, which demonstrates the user-facing functions in the ```foldcluster``` package. You can load the vignette after installation using this code.
```r
vignette("introduction", package="foldcluster")
```

## Code for Reproducing Results
* Application to the GSE81861 cell line dataset (Figures 4-6, Supplement Figure I.10): ```cells.Rmd```.
* Gaussian mixture, skew Gaussian mixture, and skew-symmetric mixture simulation studies (Supplement Figures A.2-A.4): ```SimStudy_mvn.Rmd```, ```SimStudy_mvsn.Rmd```, and ```SimStudy_nongauss.Rmd```, respectively.
* Illustration of oracle clusterings (Supplement Figures B.5 and B.6): ```convergence_plots.R```.
* Comparison of FOLD and competitors to the galaxy, acidity, enzyme, iris, flea-beetles, and wine datasets (Supplement Figures C.7-C.9): ```benchmarks.R```. Links to download the datasets are included in the comments.
* Credible ball simulation study: ```simulatespirals.R```. 

## Code for Reproducing Additional Figures
* Example of over-clustering (Figure 1): ```introexample.R```.
* Clustering the moons data (Figures 2 and 3): ```locationmoons.R```.
* Contour plots of simulation examples (Supplement Figure A.1): ```contourmaker.R```.
* Plot of the spirals data (Supplement Figure J.11): ```spiralsplot.R```.
