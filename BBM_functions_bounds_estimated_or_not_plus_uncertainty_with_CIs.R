# Set of functions to fit the Bounded Brownian Motion model
# Written by F. Boucher & V. Démery, October 7th, 2015

###############################################################################
####### SET OF AUXILIARY FUNCTIONS USED BY THE MASTER FUNCTION BELOW ##########
###############################################################################
# Diagonalize the diffusion matrix that has been discretized
# returns: the diffusion matrix, its diagonal and the base change matrix (vectors)
DiffMat=function (Npts){
	# Npts is the number of points in which the interval is discretized
	A=matrix(0,Npts,Npts)
	for (i in 1:(Npts-1)){
		A[i,i]=-2
		A[i,i+1]=1
		A[i+1,i]=1
	}
		A[1,1]=A[Npts,Npts]=-1
        eig=eigen(A)
    return(list(Diff=A,diag=diag(eig$values),passage=eig$vectors))
}  

# write to which cell of the grid a given position belongs to, 'continuous' version
VectorPos_bounds=function(x,Npts,bounds){
	X=rep(0,Npts)
	if (x==bounds[2]){X[Npts]=1}
	else {
		nx=(Npts-1)*(x-bounds[1])/(bounds[2]-bounds[1])
		ix=floor(nx)
		ux=nx-ix
		X[ix+2]=ux
		X[ix+1]=1-ux
	}	
	return(X*(Npts-1)/(bounds[2]-bounds[1]))
}

# Convolution of the matrix x with the propagator for a time t
ConvProp_bounds=function(X,t,dCoeff,dMat,bounds){
	vDiag=dMat$diag ; P=dMat$passage ; tP=t(P)
	Npts=dim(dMat$diag)[1]
	tau=((bounds[2]-bounds[1])/(Npts-1))^2
	expD=matrix(0,Npts,Npts)
	for (i in 1:Npts){expD[i,i]=exp(exp(dCoeff)*(t/tau)*diag(vDiag)[i])}
	a=P%*%expD%*%tP%*%X
	return(apply(a,1,function(x) max(x,0)))
}

# format tree and trait --> the tree is ordered from tips to root, with edge.length binded to the topology
# A list is also initiated, filled with the position (probabilistic) of each tip and 1 for internal nodes.
FormatTree_bounds=function(tree,trait,Npts,bounds){
require(ape)	
#bounds=c(min(trait),max(trait))
tree=reorder.phylo(tree,'postorder')
ntips=length(tree$tip.label)
tab=cbind(tree$edge,tree$edge.length) ; colnames(tab)=c('parent','children','brlen')
Pos=list() # one element per node
for (i in 1:(2*ntips-1)){
	if (i>ntips){
		Pos[[i]]=1
	}
	else {
		Pos[[i]]= VectorPos_bounds(trait[tree$tip.label[i]],Npts,bounds=bounds)
	}
}
return(list(tab=tab,Pos=Pos))
}

# calculate log-likelihood over the whole tree, to be minimized
LogLik_bounds=function(tree_formatted,dCoeff,dMat,bounds){
#tree_formatted obtained through FormatTree
Npts=dim(dMat$diag)[1]
tree_formatted2= tree_formatted
logFactor=0
for (i in 1:dim(tree_formatted2$tab)[1]){
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=dCoeff,dMat=dMat,bounds)
	norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
	logFactor=logFactor+log(norm)
}
return(log(max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))+logFactor)
}

# wrapper for when we only optimize dCoeff
bBM_loglik_bounds=function(tree_formatted,dMat,bounds){
	fun=function(dCoeff){
		return(-LogLik_bounds(tree_formatted,dCoeff,dMat,bounds))
	}
	return(fun)
}

##################################################
####### BOUNDS FIXED TO MIN & MAX OF TRAIT #######
##################################################
Optim_bBM_bounds_fixed=function(tree,trait,Npts=100,bounds=NULL){
	if (is.null(bounds)){
    bounds=c(min(trait),max(trait))
    }
    tree_formatted= FormatTree_bounds(tree,trait,Npts,bounds)
    dMat=DiffMat(Npts)
    fun= bBM_loglik_bounds(tree_formatted,dMat,bounds)
	opt=optim(par=var(trait)/max(branching.times(tree)),fn=fun,method='Brent',lower=-30,upper=10,hessian=FALSE)
    # dCoeff is log(sigma)
    # now retrieve the ML value of x0, using the ML of dCoeff
    tree_formatted2= tree_formatted
    for (i in 1:dim(tree_formatted2$tab)[1]){
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=opt$par,dMat=dMat,bounds)
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    }   
    x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
    res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par),root_value=x0),lnL=-opt$value,k=4,aic=2*(4+opt$value),aicc=2*(4+opt$value)+40/(length(trait)-5),method='Brent',convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	return(res)	
}

##################################################
##### NOW FUNCTIONS WITH BOUNDS OPTIMIZED TOO ####
##################################################
LogLik_bounds_opt=function(tree,trait,dCoeff,Npts,bounds){
	if ((bounds[1]>min(trait))|(bounds[2]<max(trait))) {return(-Inf)} # bounds have to be outside the trait interval
	else {
tree_formatted=FormatTree_bounds(tree,trait,Npts,bounds)
dMat=DiffMat(Npts)
tree_formatted2= tree_formatted
logFactor=0
for (i in 1:dim(tree_formatted2$tab)[1]){
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=dCoeff,dMat=dMat,bounds)
	norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
	logFactor=logFactor+log(norm)
}
return(log(max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))+logFactor)
	}
}

Optim_bBM_bounds_estimated=function(tree,trait,Npts=100){
	fun= function(pars){
		return(-LogLik_bounds_opt(tree,trait,pars[1],Npts,bounds=c(-exp(pars[2])+min(trait),exp(pars[3])+max(trait))))
	}
	opt=optim(par=c(var(trait)/max(branching.times(tree)),0,0),fn=fun,method='L-BFGS-B',lower=c(-30,-20,-20),upper=c(10,20,20),hessian=FALSE)
    # dCoeff is log(sigma)
    # now retrieve the ML value of x0, using the ML of dCoeff
    bounds=c(-exp(opt$par[2])+min(trait),exp(opt$par[3])+max(trait))
    tree_formatted=FormatTree_bounds(tree,trait,Npts,bounds)
    tree_formatted2= tree_formatted #
    dMat=DiffMat(Npts)
    for (i in 1:dim(tree_formatted2$tab)[1]){
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],dCoeff=opt$par[1],dMat=dMat,bounds)
        tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    }
    x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
    res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par[1]),root_value=x0),lnL=-opt$value,k=4,aic=2*(4+opt$value),aicc=2*(4+opt$value)+40/(length(trait)-5),method='L-BFGS-B',convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
	return(res)
	
}

##################################################
####### MASTER FUNCTION THAT FITS THE MODEL ######
##################################################
# 'tree' must be a fully dichotomous tree in phylo format
# 'trait' must be a numeric vector of trait values with names matching tree$tip.label
# 'Npts' is the number of points used for discretizing the trait interval. Default to 100.
# Smaller values should not be used, but higher (c. 200) are recommanded.
# 'bounds' specifies the bounds of the trait interval. There are three possibilities:
# 1) 'Fixed' (default) fixes bounds to the min & max of the trait vector
# Although it is not proven analytically, this is the ML estimator of the bounds
# it is also the fastest and most reliable algorithm since only 2 parameters need to be estimated: x0 and sigma^2
# 2) it is also possible to specify some values for the bounds by providing a vector with two values c(min,max). This can be useful to explore the likelihood function for the bounds around the ML estimate, but we rather recommend to use 'uncertainty=TRUE' (see below) for that purpose
# 3) If one does not trust option 1), the bounds can be estimated along with other parameters
# In our experience, this takes on average c. 10 more times than option 1) and is generally less reliable: if you use that option, always check the likelihood value returned by option 1), it should almost always be equal or higher.  
# When setting 'uncertainty=TRUE', the function will also estimate uncertainty around ML estimates of the parameter. It will calculate the confidence intervals over which the the likelihood is higher than half of the maximum likelihood and return them. 
# The 'effort_uncertainty' parameter sets the number of points that are used to calculate the likelihood surface around ML estimates (default to 100, which generally gives graphs that are smooth enough)
# In addition, 'uncertainty=TRUE' will also produce plots of:
# 1) the log-likelihood left of the lower bound, all other parameters set to their ML estimate
# 2) the log-likelihood right of the upper bound, all other parameters set to their ML estimate
# 3) the log-likelihood around the ML value of sig2, all other parameters set to their ML estimate
# 4) the density of probability for the value of the trait at the root of the tree
# in all of the 4 plots, the log-likelihood surface is plotted as a black curve and the ML value is indicated by a vertical red line

# The object returned by the function is a list containing:
# ML estimated for the 4 parameters of the BBM model ($par),
# the log-likelihood,
# the number of parameters (k=4),
# aic, aicc,
# the method used for the optimization of sigma ('Brent' or 'L-BFGS-B' depending on the option chose for 'bounds', see details in the help of the 'optim' function)
# messages from 'optim' related to possible convergence issues
# confidence intervals for parameter estimates if 'uncertainty=TRUE' (not the default)

fit_BBM_model_uncertainty=function(tree,trait,Npts=100,bounds='Fixed',uncertainty=F, effort_uncertainty=100){
	if (is.phylo(tree)==FALSE){stop("Phylogenetic tree should be in phylo format")}
	if (is.binary.tree(tree)==FALSE){stop("Phylogenetic tree should be fully dichotomous")}
	if (is.null(names(trait))){stop("Trait vector must have names")}
	if ((all(names(trait)%in%tree$tip.label)==F)|(all(tree$tip.label%in%names(trait))==F)){stop("Names of the trait vector do not match tip labels")}

if (bounds=='Fixed'){
		res=Optim_bBM_bounds_fixed(tree,trait,Npts=100,bounds=NULL)
	}
if (is.numeric(bounds)){
	if (length(bounds)==2){
		res=Optim_bBM_bounds_fixed(tree,trait,Npts=100,bounds=bounds)
		}
	else {stop("The vector of bounds should have two elements.")}	
	}
if (bounds=='Estimate'){
		res=Optim_bBM_bounds_estimated(tree,trait,Npts=100)
	}	
# uncertainty assessment
if (uncertainty==T){
	if (bounds!='Fixed'){
		stop("Uncertainty estimation possible only when bounds='Fixed'")
	}
    bounds=c(min(trait),max(trait))
    tree_formatted= FormatTree_bounds(tree,trait,Npts,bounds)
    dMat=DiffMat(Npts)
    dCoeff=log(res$par$sigsq/2) # check...
    # likelihood around lower bound
    ll_lower_bound=as.data.frame(matrix(NA, effort_uncertainty,2))
    colnames(ll_lower_bound)=c('lower_bound','loglik')
    ll_lower_bound$lower_bound=seq(from=(bounds[1]-3*(bounds[2]-bounds[1])),to=bounds[1], length.out=effort_uncertainty)
    ll_lower_bound$loglik=sapply(ll_lower_bound$lower_bound,FUN=function(x){LogLik_bounds(tree_formatted,dCoeff,dMat,bounds=c(x,bounds[2]))})
     # likelihood around upper bound
    ll_upper_bound=as.data.frame(matrix(NA, effort_uncertainty,2))
    colnames(ll_upper_bound)=c('upper_bound','loglik')
    ll_upper_bound$upper_bound=seq(from=bounds[2],to=(bounds[2]+3*(bounds[2]-bounds[1])), length.out=effort_uncertainty)
    ll_upper_bound$loglik=sapply(ll_upper_bound$upper_bound,FUN=function(x){LogLik_bounds(tree_formatted,dCoeff,dMat,bounds=c(bounds[1],x))})
    # likelihood around ML of sigma 2
    ll_sig2=as.data.frame(matrix(NA, effort_uncertainty,2))
    colnames(ll_sig2)=c('sig2','loglik')
    ll_sig2$sig2 =seq(from=(res$par$sigsq/10),to=(res$par$sigsq*10), length.out=effort_uncertainty)
    ll_sig2$loglik=sapply(ll_sig2$sig2,FUN=function(x){LogLik_bounds(tree_formatted,log(x/2),dMat,bounds)})
	# confidence intervals
	lnL=res$lnL
	CI_sigsq=c(min(ll_sig2[which(ll_sig2[,2]>(lnL-2)),1]),max(ll_sig2[which(ll_sig2[,2]>(lnL-2)),1]))
	CI_lower_bound=c(min(ll_lower_bound[which(ll_lower_bound[,2]>(lnL-2)),1]),max(ll_lower_bound[which(ll_lower_bound[,2]>(lnL-2)),1]))
	CI_upper_bound=c(min(ll_upper_bound[which(ll_upper_bound[,2]>(lnL-2)),1]),max(ll_upper_bound[which(ll_upper_bound[,2]>(lnL-2)),1]))
	cells_x0=c(min(which(res$root_density>(max(res$root_density)/2))),max(which(res$root_density>(max(res$root_density)/2))))
	CI_x0=res$par$bounds[1]+(res$par$bounds[2]-res$par$bounds[1])/(Npts-1)*(cells_x0-1)
	res$Confidence_intervals=list(CI_lower_bound = CI_lower_bound,CI_upper_bound = CI_upper_bound,CI_sigsq = CI_sigsq,CI_root_value= CI_x0)
	# plots of likelihood profile
	par(mfrow=c(2,2))
	plot(ll_lower_bound[,1:2],type='l',main='Loglik profile left of ML estimate of lower bound')
	abline(v=min(trait),col=2,lwd=2)
	plot(ll_upper_bound[,1:2],type='l',main='Loglik profile right of ML estimate of upper bound') 
	abline(v=max(trait),col=2,lwd=2)
	plot(ll_sig2[,1:2],type='l',main='Loglik profile around ML estimate of sig2',log='x') 
	abline(v=res$par$sigsq,col=2,lwd=2)	
	plot(res$root_density~seq(from=min(trait),to=max(trait),length.out=Npts),type='l',main='Density of probability for trait at the root',xlab='Root value',ylab='Probability density')
	abline(v=res$par$root_value,col=2,lwd=2)
}
	return(res)
}
