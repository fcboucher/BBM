#Below is an example of the use of BBM with simulated data:
library(ape)

# Simulate data
library(geiger) # geiger is needed for simulating the tree
tree=sim.bdtree(stop='taxa',n=20) # tree with 20 tips
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree to  a total depth of 100
TRAIT= Sim_BBM(tree,x0=0,Npts=100,sigma=1,bounds=c(-5, 5)) # TRAIT simulated on the tree, with many hits on the bounds

# Fit the model to the data simulated
BBM=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=100,bounds='Fixed',uncertainty=F) # bounds fixed to min/max of the observed trait, no assessment of uncertainty
BBM.1=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=100,bounds='Estimate',uncertainty=F) # ML value of bounds estimated along other parameters
BBM.2=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=100,bounds='Fixed',uncertainty=T, effort_uncertainty=100) # bounds fixed to min/max of the trait, but uncertainty measured 
