#COMPATIBILITY WITH TCLUST CRAN names of input and output parameters and their corresponding equivalences

#input parameters in TCLUST CRAN and equivalence here in brackets
#x A matrix or data (X)
#k The number of clusters (K)
#alpha The proportion of observations to be trimmed (alpha)
#nstart The number of random initializations to be performed 
#iter.max The maximum number of concentration steps to be performed (cstep)
#restr.fact = 1.
#equal.weights specifying whether equal cluster weights are equal (equal_weights)
#zero.tol The zero tolerance used. By default set to 1e-16 (zero.tol)


#output parameters in TCLUST CRAN and equivalence here in brackets
#centers A matrix of size p x k containing the centers columwise of each cluster. (center)
#cov An array of size p x p x k containing the covariance matrices of each cluster (sigma)
#cluster A numerical vector of size n containing the cluster assignment for each observation (assig)
#obj The value of the objective function of the best solution (obj)
#size An integer vector of size k, returning the number of observations contained by each cluster (csize)
#weights A numerical vector of length k, containing the weights of each cluster (cw)


#####FUNCTIONS DEVOTED TO APPLY CONSTRAINTS TO COVARIANCE MATRICES : restr2_eigen 
##restr2_eigen
##FUNCTION FOR APPLYING EIGEN CONSTRAINTS. These are the typical constraints
## Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26.

restr2_eigenv <- function(autovalues, ni.ini, factor_e, zero.tol)
{
  ###### function parameters:
  ###### autovalues: matrix containin eigenvalues 
  ###### ni.ini: current sample size of the clusters
  ###### factor_e: the level of the constraints
  ###### zero.tol:toletance level
  
  ###### inicializations
  c=factor_e
  d=t(autovalues)
  p = nrow (autovalues) 
  K = ncol (autovalues)	
  n=sum(ni.ini)
  nis = matrix(data=ni.ini,nrow=K,ncol=p)
  
  ###### d_ is the ordered set of values in which the restriction objective function change the definition
  ###### points in d_ correspond to  the frontiers for the intervals in which this objective function has the same definition
  ###### ed is a set with the middle points of these intervals 
  
  d_=sort(c(d,d/c))
  dim=length(d_)
  d_1=d_
  d_1[dim+1]=d_[dim]*2
  d_2=c(0,d_)
  ed=(d_1+d_2)/2
  dim=dim+1;
  
  ###### the only relevant eigenvalues are those belong to a clusters with sample size greater than 0.
  ###### eigenvalues corresponding to a clusters whit 0 individuals has no influence in the objective function
  ###### if all the eigenvalues are 0 during the smart initialization we assign to all the eigenvalues the value 1  
  
  if ((max(d[nis>0]) <= zero.tol))
    return (matrix (0, nrow = p, ncol = K))	        ##  solution corresponds to 0 matrix
  
  ###### we check if the  eigenvalues verify the restrictions
  
  if (abs(max(d[nis>0])/min(d[nis>0]))<=c)
  {
    d[nis==0]=mean(d[nis>0])
    return (t (d))					## the solution correspond to the input because it verifies the constraints
    #dfin=d
  }
  
  ###### sol[1],sol[2],.... contains the critical values of the interval functions which defines the m objective function  
  ###### we use the centers of the interval to get a definition for the function in each interval
  ###### this set with the critical values (in the array sol) contains the optimum m value  
  
  t <- s <- r <- array(0,c(K,dim))
  sol <- sal <- array(0,c(dim))
  
  for (mp_ in 1:dim) 
  {
    for (i in 1:K)
    {  
      r[i,mp_]=sum((d[i,]<ed[mp_]))+sum((d[i,]>ed[mp_]*c))
      s[i,mp_]=sum(d[i,]*(d[i,]<ed[mp_]))
      t[i,mp_]=sum(d[i,]*(d[i,]>ed[mp_]*c))
    }
    
    sol[mp_]=sum(ni.ini/n*(s[,mp_]+t[,mp_]/c))/(sum(ni.ini/n*(r[,mp_])))
    
    e = sol[mp_]*(d<sol[mp_])+d*(d>=sol[mp_])*(d<=c*sol[mp_])+(c*sol[mp_])*(d>c*sol[mp_])
    o=-1/2*nis/n*(log(e)+d/e)
    
    sal[mp_]=sum(o)
  }
  
  ###### m is the optimum value for the eigenvalues procedure 
  eo=which.max(c(sal))
  m=sol[eo]
  
  ###### based on the m value we get the restricted eigenvalues
  
  t (m*(d<m)+d*(d>=m)*(d<=c*m)+(c*m)*(d>c*m))	##	the return value
}



##      f.restr
##	function which manages constraints application 
f.restr <- function (iter, pa)
{
  
  u = array (NA, c(pa$p, pa$p, pa$K))
  d = array (NA, c(pa$p, pa$K))
  
  for (k in 1:pa$K)
  {
    ev = eigen (iter$sigma[,,k])
    u [,,k] <- ev$vectors
    d [,k] <- ev$values
  }
  
  d [d < 0] <- 0		##	all eigenvalue < 0 are assigned to 0, this issue appears for numerical errors
  
  d=restr2_eigenv (autovalues=d, ni.ini=iter$csize, factor_e=pa$restr.fact, zero.tol=pa$zero.tol) 
  ##	checking for singularity in all clusters.
  iter$code = max(d) > pa$zero.tol
  
  if (!iter$code) return (iter)
  
  for (k in 1:pa$K)	##	reconstructing the sigma matrices
    iter$sigma[,,k] <- u[,,k] %*% diag (d[,k], nrow = pa$p) %*% t(u[,,k])
  
  return (iter)
}

###### MISCELANEOUS FUNCTIONS: dmnorm ssclmat TreatSingularity

##	Multivariate normal density
dmnorm <- function(X,mu,sigma) ((2*pi)^(-length(mu)/2))*(det(sigma)^(-1/2))*exp(-0.5*mahalanobis(X,mu,sigma))

##	get a matrix object out of the sigma 
ssclmat <- function (x, k) as.matrix (x[,,k])

##	to manage simgular situations 
TreatSingularity <- function (iter, pa) 
{	
  warning ("Data in no general position")
  return (iter)
}


######## FUNCTIONS FOR RANDOM STARTS:  getini    InitClusters
##	calculates the initial cluster sizes
getini <- function (K, no.trim)
{
  if (K == 1)
    return (no.trim)
  
  pi.ini  <- runif(K)
  ni.ini <- sample(x = K, size = no.trim, replace = T, prob = pi.ini / sum (pi.ini))
  return (tabulate(ni.ini, nbins = K))
}

##	calculates the initial cluster assignment and initial values for the parameters
InitClusters <- function (X, iter, pa)
{
  dMaxVar = 0
  for (k in 1:pa$K)
  {
    idx <- sample (1:pa$n, pa$p+1)
    X.ini = X [drop = F,idx,]#sample (1:pa$n, pa$p+1),]	##	selecting observations randomly for the current init - cluster
    
    iter $center [k,] <- colMeans (X.ini)			##	calculating the center
    cc <- (pa$p/(pa$p+1))*cov (X.ini)			##	calculating sigma (cc = current cov)
    iter$sigma[,,k] <- cc
  }
  
  if (pa$equal.weights)						##	if we're considering equal weights, cw is set here AND NEVER CHANGED
  { iter$csize <- rep (pa$no.trim / pa$K, pa$K)}	else	{iter$csize = getini (pa$K, pa$no.trim)}
  iter$cw <- iter$csize / pa$no.trim				## if we're considering different weights, calculate them, and they're gonna be recalculated every time in findClustAss
  return (iter)
}


######## FUNCTION FOR estimating model parameters: estimClustPar
##PAPERS
##HARD ASSIGNMEMT Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26. 
##MIXTURE MODEL Garcia-Escudero, Luis Angel, Alfonso Gordaliza, and Agustin Mayo-Iscar. "A constrained robust proposal for mixture modeling avoiding spurious solutions." Advances in Data Analysis and Classification 8.1 (2014): 27-43.

estimClustPar <- function (X, iter, pa)
{				
  
  for (k in 1:pa$K)
  {
    if (iter$csize[k] > pa$zero.tol)	##	if cluster's size is > 0
    {
      iter$center[k,] = (t(iter$z_ij[,k]) %*% X) / iter$csize[k]
      X.c <- (X - matrix (iter$center[k,], ncol = pa$p, nrow = pa$n, byrow = TRUE))
      iter$sigma[,,k] <- (t(X.c * iter$z_ij[,k]) %*% X.c) / iter$csize[k]
    }			else		iter$sigma[,,k] <- 0	
  }
  
  return (iter)
}


######## FUNCTIONS FOR obtaining the assigment and trimming: findClustAssig (mixture models and hard assigment) 
##PAPERS
##HARD ASSIGNMEMT Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26. 
##MIXTURE MODEL Garc?a-Escudero, Luis Angel, Alfonso Gordaliza, and Agust?n Mayo-Iscar. "A constrained robust proposal for mixture modeling avoiding spurious solutions." Advances in Data Analysis and Classification 8.1 (2014): 27-43.
findClustAssig <- function (X, iter, pa)             
{										
  ll = matrix (NA, pa$n, pa$K)	
  
  for (k in 1:pa$K)	ll[,k] <- iter$cw[k] * dmnorm(X,iter$center[k,],ssclmat (iter$sigma,k)) ## dmvnorm could be used here...
  
  old.assig <- iter$assig
  iter$assig <- apply(ll,1,which.max)						##	searching the cluster which fits best for each observation
  
  pre.z_h=apply(ll,1,'max')
  
  pre.z_m=apply(ll,1,'sum')
  
  pre.z_=matrix(pre.z_m, nrow=pa$n, ncol=pa$K,byrow=FALSE)
  
  
  ##### To obtain the trimming:  tc.set is the non trimming indicator
  if (pa$opt=="MIXTURE") tc.set=(rank(pre.z_m, ties.method="random")> floor(pa$n*(pa$alpha)))
  if (pa$opt=="HARD")    tc.set=(rank(pre.z_h, ties.method="random")> floor(pa$n*(pa$alpha)))
  
  #hard assigment including trimming
  iter$assig <- apply(ll,1,which.max)*tc.set
  
  ##### To obtain the iter$z_ij matrix contining the assigment and trimming
  if (pa$opt=="MIXTURE") iter$z_ij=ll/(pre.z_ + (pre.z_ ==0))*tc.set   #mixture assigment including trimming 
  if (pa$opt=="HARD")
  {
    iter$z_ij=0*iter$z_ij
    iter$z_ij[cbind ((1:pa$n), iter$assig+(iter$assig==0)  )] <- 1	
    iter$z_ij[tc.set==FALSE,]=0
  }
  
  
  if (pa$opt=="HARD")      iter$code <- 2 * all (old.assig == iter$assig)		##	setting the code - parameter, signaling whether the assignment is the same than the previous --- is the only stopping rule implemented
  
  ##### To obtain the size of the clusters and the estimated weight of each population
  if (pa$opt=="HARD")	  iter$csize <- tabulate (iter$assig, pa$K)
  if (pa$opt=="MIXTURE")    iter$csize <- apply(iter$z_ij, 2,sum)  ##/(1-pa$alpha) 
  
  if (!pa$equal.weights)     iter$cw <- iter$csize / sum(iter$csize)					##	and cluster weights
  
  
  return (iter)
}


######## FUNCTION FOR obtaining the objetive functions value for mixture (obj_m) hard (obj_h) 
##HARD ASSIGNMEMT Fritz, H., Garcia-Escudero, L. A., & Mayo-Iscar, A. (2012). tclust: An r package for a trimming approach to cluster analysis. Journal of Statistical Software, 47(12), 1-26. 
##MIXTURE MODEL Garc?a-Escudero, Luis Angel, Alfonso Gordaliza, and Agust?n Mayo-Iscar. "A constrained robust proposal for mixture modeling avoiding spurious solutions." Advances in Data Analysis and Classification 8.1 (2014): 27-43.

calcobj <- function (X, iter, pa) 
{			
  ww_m=matrix(0,nrow=pa$n,ncol=1)
  
  ww_h=matrix(0,nrow=pa$n,ncol=1)
  
  
  for (k in 1:pa$K) {
    
    
    w_m=iter$cw[k]*dmnorm(X,iter$center[k,],ssclmat (iter$sigma,k)) 
    
    
    ww_m=w_m*(w_m>=0)+ww_m     ##	calculates each individual contribution for the obj funct mixture
    
    w_h=w_m*(iter$assig==k)
    
    ww_h=w_h*(w_h>=0)+ww_h     ##	calculates each individual contribution for the obj funct hard
    
  }
  
  ww_m=ww_m*(ww_m>=0)
  
  ww_h=ww_h*(ww_h>=0)
  
  
  if (pa$opt=="MIXTURE")		        iter$obj <- sum(log(ww_m[iter$assig>0]))  
  
  if (pa$opt=="HARD")			iter$obj <- sum(log(ww_h[iter$assig>0]))  
  
  return (iter)
  
}

######## AUXILIAR FUNCTION TO perfom the concentration steps max.iter times
tclust_update <- function(X, pa, iter, max.iter){
  
  for (i in 0:max.iter)
  {
    ##	restricting the clusters' scatter structure
    iter <- f.restr(iter=iter,pa=pa)	
    
    ##  estimates the cluster's assigment and TRIMMING (MIXT models and HARD)
    iter <- findClustAssig(X,iter,pa)	 	
    
    if ((iter$code==2) ||		##	if findClustAssig returned code=2 (convergence)
        (i == max.iter))		##	or we're in the last concentration step:
      break			##	break the for - loop - we finished this iteration!
    ##  don't re-estimate cluster parameters this time
    
    iter <- estimClustPar(X, iter, pa)		## estimates the cluster's parameters 
  }
  
  iter <- calcobj (X, iter, pa)			##  calculates the objective function value
  
  return(iter)
}



######## MAIN FUNCTION ####################################

# install.packages("foreach")
# install.packages("iterators")
library(doRNG)
library(doParallel)

## List of auxiliary functions to be exported to parallel workers
f_aux_tclust <- c( "restr2_eigenv","f.restr", "getini", "InitClusters", 
                   "estimClustPar", "findClustAssig", "calcobj", "dmnorm", "ssclmat", 
                   "tclust_update","TreatSingularity")

## tclust Parallelized Main Function

tclust <- function(X,K,alpha = 0.05, nstart = 20, cstep1=2, cstep2=10, nkeep=5,  
                   equal.weights = FALSE, restr.fact=5,zero.tol = 1e-16,  trace = 0,   opt="HARD",
                   sol_ini_p = FALSE, sol_ini = NA )  {
  
  # --- Data Validation ---
  if (!is.numeric (X)) 
    stop ("parameter x: numeric matrix/vector expected")
  if( !is.matrix (X))
    X <- matrix (X, ncol = 1)
  
  n <- nrow (X)
  p <- ncol (X)
  no.trim <- floor(n*(1-alpha))
  
  # preparing lot's of lists: pa (input parameters) and iter (current value of the parameters)
  
  # --- Input Parameters List (pa) ---
  pa <- list (		            ##	 input parameters for the procedure
    n = n,			    ##	number of observations
    p = p,	                    ##	number of dimensions
    alpha=alpha,		    ##	level of trimming
    trimm = n-no.trim,	    ##	number of observations which are considered as to be outliers
    no.trim = no.trim,	    ##	number of observations which are considered as to be not outliers
    K = K,			    ##	number of clusters 
    equal.weights = equal.weights,		##	equal population proportions  for all clusters
    zero.tol = zero.tol,			##	zero tolerance	(to be implemented)	
    restr.fact=restr.fact,                    ##	level eigen constraints
    opt = opt,                              ##	estimated model "MIXTURE" (mixture model) "HARD" (hard assignment)  "FUZZY" (fuzzy)	
    sol_ini_p = sol_ini_p,                  ##	initial solution provided by the user TRUE/FALSE   if TRUE is stored in sol_ini
    cstep1=cstep1,                          ##	number of iterations
    cstep2=cstep2,
    nstart=nstart,
    nkeep=nkeep
  )
  
  # --- Internal State List (iter) ---
  iter <- list (					##	current value of the parameters
    obj = NA,				##	current objective value
    assig = array (0, n),			##	cluster group of assignment
    csize = array (NA, K),			##	cluster number of observations assigned to the cluster 
    cw = rep (NA, K),			##	cluster estimated proportions 
    sigma = array (NA, c (p, p, K)),	##	cluster's covariance matrix
    center = array (NA, c(K, p)),		##	cluster's centers
    code = NA,				##	this is a return code supplied by functions like findClustAssig
    z_ij = matrix (0, nrow=n, ncol=K )	##	cluster assignment given by 0/1 columns               
  )
  
  # --- Parallel Backend Initialization ---
  # Adjust threads for OpenMP if necessary and setup the cluster
  Sys.setenv(OMP_NUM_THREADS="4")
  parclus <- makeCluster(20)
  registerDoParallel(parclus)
  on.exit(stopCluster(parclus))
  
  
  ##############################################################################
  ######### PHASE 1: INITIALIZATIONS AND SHORT CONCENTRATION STEPS #############
  ##############################################################################
  
  if (pa$sol_ini_p==TRUE) 
  {
    ##### Single-start execution if an initial solution is provided
    pa$nstart=1 ; iter=sol_ini
    return(tHDDC_update(X, pa, iter, max.iter=pa$cstep1))
  }
  else{
    ##### Generate multiple random starts and run for cstep1 iterations
    iniciales <- foreach(inicial=1:pa$nstart, .export=f_aux_tclust, .inorder=TRUE) %dorng% {
      return(tclust_update(X, pa, InitClusters(X,iter,pa),pa$cstep1))}
  }
  
  
  print("j")  
  
  ##############################################################################
  ########## PHASE 2: SELECT THE nkeep BEST SOLUTIONS (BY BIC)  ################
  ##############################################################################
  
  ##### Sort solutions based on BIC (lower is better) and pick top candidates
  mejores <- order(sapply(iniciales, function(iter) iter$obj), decreasing=TRUE)[1:pa$nkeep]  
  
  print(mejores)
  
  ##############################################################################
  ########## PHASE 3: ITERATE TOP SOLUTIONS UNTIL CONVERGENCE  #################
  ##############################################################################
  
  ##### Run the selected best candidates for a longer duration (cstep2)
  finales <- foreach(l=mejores, .export=f_aux_tclust) %dorng% {
    return(tclust_update(X, pa, iniciales[[l]], pa$cstep2))} 
  
  ##### Find the overall best solution from the refined results
  j <- which.max(sapply(finales, function(iter) iter$obj))
  
  ##### Final output
  return(finales[[j]])
  
}
