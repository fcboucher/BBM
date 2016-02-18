#Below is an example of the use of BBM with simulated data:
library(ape)

# Simulate data
library(geiger) # geiger is needed for simulating the tree
tree=sim.bdtree(stop='taxa',n=50) # tree with 50 tips
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree to a total depth of 100
TRAIT= Sim_BBM(tree,x0=0,Npts=50,sigma=1,bounds=c(-5, 5)) # TRAIT simulated on the tree, with many hits on the bounds: for that you need to source the function 'Sim_BBM.R'
hist(TRAIT,breaks=20) # the distribution of the trait at the tips of the tree: it should be rather flat...

# Fit the model to the data simulated: we use only 50 points for discretizing the trait interval to make it faster, but more points should be used on empirical datasets
BBM=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=50,bounds='Fixed',uncertainty=F) # bounds fixed to min/max of the observed trait, no assessment of uncertainty
BBM.1=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=50,bounds='Estimate',uncertainty=F) # ML value of bounds estimated along other parameters
BBM.2=fit_BBM_model_uncertainty(tree,trait=TRAIT,Npts=50,bounds='Fixed',uncertainty=T, effort_uncertainty=100) # bounds fixed to min/max of the trait, but uncertainty measured. This will produce a plot of log-likelihood profiles around ML estimates of parameters 

BBM$par # parameters estimated
BBM.1$par # in particular, check the ML estimates of the bounds: are they close to the ones fixed in the previous run (BBM)?

# which optimization performed best?
BBM$lnL
BBM.1$lnL # normally, BBM.1 should have a likelihood higher or equal than BBM, since BBM is a special case of BBM.1 (we fixed bounds in BBM). However, it often happens that BBM gets a higher likelihood, this is because it is easier to optimize since it has fewer parameters...

# AICc comparison with other classic models of evolution implemented in package {geiger}
BM=fitContinuous(phy=tree,dat=TRAIT,model='BM') # Brownian motion with no bounds
OU=fitContinuous(phy=tree,dat=TRAIT,model='OU') # Ornstein-Uhlenbeck process with a single optimum

BBM$aicc # does BBM get the lowest AICc?
BM$opt$aicc
OU$opt$aicc

# output the characteristic time for the trait to cross the trait interval (i.e. evolve from one bound to the other)
tau=(BBM$par$bounds[2]-BBM$par$bounds[1])^2/BBM$par$sigsq
tau # this is the characteristic time: it roughly compares the width of the trait interval (numerator) to the rate of diffusion (denominator). See our paper for more details.
max(branching.times(tree))/tau # total tree depth divided by the characteristic time of the process. This is the expected number of crossings of the trait interval (caution: not the number of hits on the bounds, many of these can happen in a row if the trait stays close to one bound). If this measure is higher than 1 then trait evolution has been strongly bounded. See our paper for more details.


