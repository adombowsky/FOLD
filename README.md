# FOLD
This is the GitHub repository for [Bayesian Clustering via Fusing of Localized Densities](https://arxiv.org/abs/2304.00074). In this article, we present Fusing of Localized Densities (FOLD), a decision theoretic clustering method that focuses on grouping the _localized densities_ from a mixture model. Here, you will find the code required to reproduce all plots, the simulation studies, and the application to the cell line dataset. Here is a quick summary of what files are included in this repository.

* ```SimStudy_mvn.Rmd```, ```SimStudy_mvsn.Rmd```, and ```Simstudy_nongauss.Rmd```: simulation studies for the Gaussian mixture, skew Gaussian mixture, and skew-symmetric mixture.
* ```simulatespirals.R```: the simulation study focused on credible balls from the Supplementary.
* ```cells.Rmd```: the cell line dataset application. [Click here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81861) for more information on the data.
* ```locationmoons.R```: an in-depth example focused on fitting a location mixture of Gaussians to the crescent moons data. 
* Other files correspond to figures used in the main text and Supplementary.
