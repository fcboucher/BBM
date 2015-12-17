#Below is an example of the use of BBM with simulated data:
library(ape)

# Simulate data
library(geiger) # geiger is needed for simulating the tree
tree=sim.bdtree(stop='taxa',n=50) # tree with 20 tips
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree to  a total depth of 100
TRAIT= Sim_BBM(tree,x0=0,Npts=50,sigma=1,bounds=c(-5, 5)) # TRAIT simulated on the tree, with many hits on the bounds: for that you need to source the function 'Sim_BBM.R'
hist(TRAIT,breaks=20) # the distribution of the trait at the tips of the tree is rather flat...

# Fit the model to the data simulated
BBM=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=50,bounds='Fixed',uncertainty=F) # bounds fixed to min/max of the observed trait, no assessment of uncertainty
BBM.1=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=50,bounds='Estimate',uncertainty=F) # ML value of bounds estimated along other parameters
BBM.2=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=50,bounds='Fixed',uncertainty=T, effort_uncertainty=100) # bounds fixed to min/max of the trait, but uncertainty measured. This will produce a plot of log-likelihood profiles around ML estimates of parameters 

BBM$par # parameters estimated
BBM.1$par # check the ML estimates of the bounds

# Comparison with other classic models of evolution implemented in package {geiger}
BM=fitContinuous(phy=tree,dat=TRAIT,model='BM') # Brownian motion with no bounds
OU=fitContinuous(phy=tree,dat=TRAIT,model='OU') # Ornstein-Uhlenbeck process with a single optimum

BBM$aicc
BM$opt$aicc
OU$opt$aicc # the model with the lowest AICc is the 'best'
