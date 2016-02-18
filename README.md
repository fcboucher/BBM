BBM provides a set of R functions to fit a bounded brownian motion (BBM) model of evolution to phylogenetic comparative data using maximum likelihood. In the BBM model, a continuous trait evolves under brownian motion with a constant evolutionary rate, but between two reflective bounds.

<img width="70%" alt="capture d ecran 2016-02-18 a 22 44 39" src="https://cloud.githubusercontent.com/assets/15023761/13159291/4450bcdc-d691-11e5-8138-bcb94d22a0df.png">

These functions only depend on the {ape} package in R and likelihoods are compatible with those of other models fitted by the 'fitContinuous' function in package {geiger}.

Functions were written by [Florian Boucher](https://sites.google.com/site/floriaboucher/) and [Vincent Démery](https://www.pct.espci.fr/~vdemery/). You can find details on the maths and code in the manuscript 'Boucher&Démery_Main_text_plus_Appendix.pdf' (accepted in Systematic Biology). 

In order to use the code, copy or download the 'BBM_functions_bounds_estimated_or_not_plus_uncertainty_with_CIs.R' R script, send all of it to R, and finally use the master function: 'fit_BBM_model_uncertainty'. Details on its use are given as comments in the script: please read them.

An example of use is given in the script 'Example.R'. It uses the 'Sim_BBM.R' function to simulate traits evolving under BBM.
