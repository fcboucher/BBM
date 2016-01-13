BBM provides a set of R functions to fit a bounded brownian motion (BBM) model of evolution to phylogenetic comparative data using maximum likelihood. In the BBM model, a continuous trait evolves under brownian motion with a constant evolutionary rate, but between two reflective bounds.

These functions only depend on the {ape} package in R and likelihoods are compatible with those of other models fitted by the 'fitContinuous' function in package {geiger}.

Functions were written by Florian Boucher and Vincent Démery. You can find details on the maths and code in the manuscript 'Boucher&Démery_Main_text_plus_Appendix.pdf' (accepted pending minor revisions in Systematic Biology). 

In order to use the code, copy or download the 'BBM_functions_bounds_estimated_or_not_plus_uncertainty_with_CIs.R' R script, send all of it to R, and finally use the master function: 'fit_BBM_model_uncertainty'. Details on its use are given as comments in the script: please read them.

An example of use is given in the script 'Example.R'. It uses the 'Sim_BBM.R' function to simulate traits evolving under BBM.
