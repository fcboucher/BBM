BBM provides a set of R functions to fit a bounded brownian motion (BBM) model of evolution to phylogenetic comparative data using maximum likelihood. In the BBM model, a continuous trait evolves under brownian motion with a constant evolutionary rate, but between two reflective bounds.

These functions only depend on the {ape} package in R and likelihoods are compatible with those of other models fitted by the 'fitContinuous' function in package {geiger}.

Functions were written by [Florian Boucher](https://sites.google.com/site/floriaboucher/) and [Vincent Démery](https://www.pct.espci.fr/~vdemery/). You can find details on the maths and code in the [preprint manuscript](https://github.com/fcboucher/BBM/blob/master/Boucher%26De%CC%81mery_Main_text_plus_Appendix.pdf) or read the published version of the article presenting the method:
*Inferring bounded evolution in phenotypic characters from phylogenetic comparative data, F.C. Boucher & V. Démery, 2016, Systematic biology, 65(4):651-661.*

In order to use the code, copy or download the 'BBM_functions_bounds_estimated_or_not_plus_uncertainty_with_CIs.R' R script, send all of it to R, and finally use the master function: 'fit_BBM_model_uncertainty'. Details on its use are given as comments in the script: please read them.

There are currently no help pages for these functions but R scripts are heavily commented and should give you all the information needed on parameters, outputs, etc.

An example of use is given [here](https://github.com/fcboucher/BBM/blob/master/R/Example.R). It uses the 'Sim_BBM.R' function to simulate traits evolving under BBM.
