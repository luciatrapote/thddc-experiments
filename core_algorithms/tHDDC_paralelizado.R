#### CATELL'S METHOD FUNCTION

catell <- function(d,thres)
{
  ###### function parameters:
  ###### d: vector containing ordered eigenvalues
  ###### thres: threshold for Catell's scree test
  dif <- abs(diff(d))
  return(max(which(dif > thres)))
}


#####FUNCTIONS DEVOTED TO APPLY CONSTRAINTS TO COVARIANCE MATRICES : restr2_eigen 

restr2_eigenv_tHDDC  <- function(autovalues, ni.ini, factor_e, zero.tol)
{
  ###### function parameters:
  ###### autovalues: matrix containing eigenvalues
  ###### ni.ini: current sample size of the clusters
  ###### factor_e: the level of the constraints
  ###### zero.tol:tolerance level
  
  ###### Inicializations
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
  
  ###### The only relevant eigenvalues are those belong to a clusters with sample size greater than 0.
  ###### eigenvalues corresponding to a clusters whit 0 individuals has no influence in the objective function
  ###### if all the eigenvalues are 0 during the smart initialization we assign to all the eigenvalues the value 1
  if ((max(d[nis>0]) <= zero.tol))
    return (matrix (0, nrow = p, ncol = K))	        ## Returns zero matrix solution
  
  ###### Check if the  eigenvalues verify the restrictions
  if (abs(max(d[nis>0])/min(d[nis>0]))<=c)
  {
    d[nis==0]=mean(d[nis>0])
    return (t (d))					## Returns input as it already satisfies constraints
  }
  
  ###### Build the 'sol' array
  ###### sol[1],sol[2],.... this array contains the critical values of the interval functions which defines the m objective function
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
  
  ###### 'm' is the optimum value for the eigenvalues procedure
  eo=which.max(c(sal))
  m=sol[eo]
  
  ###### Return restricted eigenvalues based on the optimal 'm'
  t (m*(d<m)+d*(d>=m)*(d<=c*m)+(c*m)*(d>c*m))	##	the return value
}




# Function to adjust the eigenvalues by groups (ensures order is maintained after restriction)

adjust_eigenvalues <- function(matriz, small, q, p) {
  ###### function parameters:
  ###### matriz: matrix containing the eigenvalues of the intrinsic dimensions
  ###### small: vector containing the noise eigenvalues of residual spaces
  ###### q: vector containing the current intrinsic dimensions
  ###### p:total dimension of the space
  
  cambios <- FALSE  # Helper variable to track if changes occurred
  
  contador <-0
  
  for (j in which(q!=0)) {
    
    contador <- contador + 1
    t <- which(matriz[1:q[j], contador] < small[contador])[1]
    
    if (!is.na(t)&& is.na(t)<=q[j]) {
      # If an intrinsic eigenvalue is smaller than the noise level, they are averaged
      valores <- c(matriz[t:q[j], contador], rep(small[contador], p-q[j]))
      nuevo_valor <- mean(valores)
      
      matriz[t:q[j], contador] <- nuevo_valor
      small[contador] <- nuevo_valor
      
      cambios <- TRUE 
    } 
  }
  
  list(matriz = matriz, small = small, cambios = cambios)
}

##	Function managing the application of constraints

f.restr_tHDDC <- function (iter, pa)
{
  ###### Replace negative eigenvalues with a near-zero floor
  iter$evalues [iter$evalues < 0] <- 1e-16
  
  repetir <- TRUE
  max_iter <- 10  # Safety limit for the adjustment loop
  
  ###### Matrix to store large (intrinsic) eigenvalues, excluding clusters with q=0
  matriz<-matrix(0,pa$p,sum(iter$q != 0))
  cont<-0
  for (i in which(iter$q != 0)){
    cont<-cont+1
    matriz[1:iter$q[i],cont]<-iter$evalues[1:iter$q[i],i]
  }
  
  ###### Vector to store small (noise) eigenvalues, excluding empty residual spaces
  small<-rep(0,sum( (pa$p-iter$q) != 0))
  cont<-0
  for (i in which( (pa$p-iter$q) != 0)){
    cont<-cont+1
    small[cont]<-iter$evalues[pa$p,i]
  }
  
  cont <- 0
  ###### Iterative process to enforce constraints across both intrinsic and residual spaces
  while (repetir && cont < max_iter) {
    cont <- cont + 1
    ###### Apply constraints on eigenvalues (large vs small)
    aux <- restr2_eigenv_tHDDC(cbind(matriz[matriz>0]), ni.ini = rep(1,sum(iter$q)), factor_e = pa$restr.fact[1], zero.tol = 1e-16)
    small <- restr2_eigenv_tHDDC(cbind(small), ni.ini = pa$p-iter$q, factor_e = pa$restr.fact[2], zero.tol = 1e-16)
    
    ###### Reconstruct the adjusted matrix
    cont_aux2 <- 1
    cont_aux <- 0
    for (i in which(iter$q != 0)){
      cont_aux <- cont_aux + 1
      matriz[1:iter$q[i], cont_aux]<-aux[cont_aux2:(cont_aux2+iter$q[i]-1)]
      cont_aux2 <- cont_aux2 + iter$q[i]
    }
    
    ###### Order constraint adjustment within groups
    resultado <- adjust_eigenvalues(matriz, small, iter$q, pa$p)
    
    matriz <- resultado$matriz
    small <- resultado$small
    
    ###### If there are changes, set repeat=TRUE and continue iterating through the loop
    repetir <- resultado$cambios  
    
  }
  
  ###### Map results back to the iter$evalues structure
  cont<-0
  for (i in which(iter$q != 0)){
    cont<-cont+1
    iter$evalues[1:iter$q[i],i]<-matriz[1:iter$q[i],cont]
  }
  
  cont<-0
  for (i in which( (pa$p-iter$q) != 0)){
    cont<-cont+1
    iter$evalues[(iter$q[i]+1):pa$p,i]<-rep(small[cont],pa$p-iter$q[i])
  }
  
  ######	Check for numerical singularity; code=1 if valid variability exists
  iter$code <- 1*( sum(small) > pa$zero.tol)
  if ((!iter$code)) {
    print("zero variability in the spherical dimensions")

  }
  return (iter)
}


######## FUNCTIONS FOR RANDOM STARTS:  getini_tHDDC    InitClusters_tHDDC

##	Calculates the initial cluster sizes (number of observations per group)

getini_tHDDC <- function (K, no.trim)
{
  if (K == 1) return(no.trim)
  
  ###### Generate random mixing proportions
  pi.ini  <- runif(K)
  
  ###### Sample initial assignments based on these proportions
  ni.ini <- sample(x = K, size = no.trim, replace = T, prob = pi.ini / sum (pi.ini))
  return(tabulate(ni.ini, nbins = K))
}


##	Calculates the initial cluster assignment and initial parameter values

InitClusters_tHDDC <- function (X, iter, pa)
{
  d <- matrix(0,nrow=pa$p,ncol=pa$K)
  v <-list()
  d0 <- rep(0,pa$K)
  
  for (k in 1:pa$K)
  {
    ###### Randomly sample q_k + 2 points to avoid singularity during initialization
    idx <- sample (1:pa$n, pa$q[k]+2)
    X.ini = X [drop = F,idx,]	
    
    ###### INITIAL CENTERS (Means)
    iter$center[k,] <- colMeans (X.ini)		
    
    ###### INITIAL COVARIANCE STRUCTURE 
    X.ini.cen <- X.ini-rep(1,pa$q[k]+2) %*% t(iter$center[k,]) 
    
    matrix <- (1/(pa$q[k]+1)) * X.ini.cen %*% t(X.ini.cen)  
    
    if (pa$q[k]>0){
      ev <- eigen(matrix)
      v[[k]] <- matrix(ev$vectors[,1:pa$q[k]], ncol = pa$q[k]) 
      d[1:pa$q[k],k] <- ev$values[1:pa$q[k]] 
      u<-array(0, dim = c(pa$p,pa$q[k],pa$K))
      for (j in 1:pa$q[k]){
        u[,j,k] <- (1/sqrt(pa$q[k]+1))*(1/sqrt(d[j,k]))*t(X.ini.cen)%*%v[[k]][,j] 
      }
      iter$evectors[[k]] <-u[,,k]
      ###### Estimate noise variance (residual subspace)
      d0[k] <- (sum(diag(cov(X.ini)))-sum(d[1:pa$q[k],k]))/(pa$p - pa$q[k])
    }
    else
    {
      ###### Case where intrinsic dimension q = 0
      d0[k] <- sum(diag(matrix))/pa$p
    }
  }
  
  ###### Fill the residual eigenvalues
  for (k in 1:pa$K){
    d[(pa$q[k]+1):pa$p,k]<-rep(d0[k],length((pa$q[k]+1):pa$p))
  }
  
  iter$evalues <- d
  
  ###### INITIAL MIXING WEIGHTS (PROPORTIONS)
  ###### If we're considering equal weights, cw is set here AND NEVER CHANGED
  if (pa$equal.weights)  {iter$csize <- rep (pa$no.trim/pa$K, pa$K)}  else {iter$csize = getini_tHDDC(pa$K, pa$no.trim)}
  
  iter$cw <- iter$csize/pa$no.trim
  
  iter$q=pa$q
  
  return (iter)
}


######## FUNCTION FOR estimating model parameters: estimClustPar_tHDDC

estimClustPar_tHDDC <- function (X, iter, pa)
{				
  
  d <- matrix(NA, pa$p, pa$K)
  aux_evectors <- list()
  
  ###### Determine whether q is fixed or estimated via Catell's method
  if (pa$catell_Y==TRUE) iter$q=array(0,pa$K)  else iter$q=pa$q
  
  iter$evalues=array(NA,c(pa$p,pa$K))
  for (k in 1:pa$K)
  {
    if (iter$csize[k] > pa$zero.tol)	##	if cluster's size is > 0
    {
      ###### Update cluster center (Weighted mean)
      iter$center[k,] <- (t(iter$z_ij[,k]) %*% X) / iter$csize[k]
      
      ###### Update weighted covariance matrix
      X.c <- (X - matrix(iter$center[k,], ncol = pa$p, nrow = pa$n, byrow = TRUE))
      cov.wei <- (t(X.c * iter$z_ij[,k]) %*% X.c) / iter$csize[k]
      
      ###### Determine if eigenvalue calculation is necessary
      if (pa$catell_Y==TRUE) calc_eig=TRUE else {if (iter$q[k]>0) calc_eig=TRUE else calc_eig=FALSE} 
      
      ###### Compute eigenvalues/eigenvectors
      if (calc_eig==TRUE) {
        ev <- eigen(cov.wei)
        d[,k] <- ev$values
        d [d[,k] < 0,k] <- 0	}
      
      ###### Automatic intrinsic dimension selection via Catell's method
      if (pa$catell_Y==TRUE) iter$q[k]= max(min(catell(d[,k],pa$catell_threshold),floor(pa$dim_max)),0)
      
      if (iter$q[k]>0){
        iter$evectors[[k]] <- ev$vectors[,1:iter$q[k]]
        iter$evalues[1:iter$q[k],k] <-d[1:iter$q[k],k]
        ###### Average variance for the residual subspace (spherical part)
        iter$evalues[(iter$q[k]+1):pa$p,k] <-(1/(pa$p-iter$q[k]))*(sum(diag(cov.wei))-sum(d[1:iter$q[k],k]))}
      else {
        ###### Purely spherical case
        iter$evalues[(iter$q[k]+1):pa$p,k] <-(1/pa$p)*sum(diag(cov.wei))}
      
    }	 else{	
      ###### Handle empty clusters by resetting to default variability
      iter$evalues[,k] <-rep(1,pa$p)}
  }
  
  return (iter)
}


####### OBJECTIVE FUNCTION COMPONENTS: A.fun
# Computes distances in both intrinsic and residual subspaces

A.fun <- function(X,pa,iter,k)
{
  if(iter$q[k]>0)
  {
    A <- matrix(iter$evectors[[k]],ncol=iter$q[k])
    Xcen <- X - rep(1,pa$n) %*% t(iter$center[k,])
    TT <- Xcen %*% A
    ###### Squared Mahalanobis distance in intrinsic space (A1)
    A1 <- rowSums( (TT^2) %*% diag((iter$evalues[1:iter$q[k],k])^(-1),iter$q[k],iter$q[k]) )
    ###### Squared distance to the intrinsic space (Orthogonal Distance) (A2)
    A2 <- ( ( iter$evalues[iter$q[k]+1,k] )^(-1) ) * rowSums( (Xcen %*% (diag(pa$p) - A %*% t(A)))^2)
    ###### Log-determinant contribution of intrinsic dimensions (A3)
    A3 <- sum( log(iter$evalues[1:iter$q[k],k]) )
  }
  else
  { 
    A1 <- A3 <- 0 
    Xcen <- X - rep(1,pa$n) %*% t(iter$center[k,])
    ###### Purely spherical distance (A2)
    A2<-(iter$evalues[iter$q[k]+1,k] )^(-1) * apply(Xcen, 1, function(x) sum(x^2))
  }
  ###### Log-determinant contribution of the residual space (A4)
  A4 <- ( pa$p - iter$q[k] ) * log( iter$evalues[iter$q[k]+1,k] )
  ###### Constant term (A5)
  A5 <- pa$p * log(2*pi)
  return(list(A1=A1, A2=A2, A3=A3, A4=A4, A5=A5))
}

######## FUNCTION FOR obtaining the assigment and trimming: findClustAssig_tHDDC

findClustAssig_tHDDC <- function(X,iter,pa)
{	
  ll <- matrix (NA, pa$n, pa$K)
  for (k in 1:pa$K){
    A <- A.fun(X, pa, iter, k)
    ll[,k] <- iter$cw[k]*exp(- 0.5 * ( A$A1 + A$A2 + A$A3 + A$A4 + A$A5 )) 
  }
  old.assig <- iter$assig
  
  pre.z_h <- apply(ll,1,max) # if we are working with hard assignament
  pre.z_m <- apply(ll,1,sum) # if we are working with MIXT assignament
  pre.z_ <- matrix(pre.z_m, nrow=pa$n, ncol=pa$K,byrow=FALSE)
  
  ##### TRIMMING STEP: tc.set identifies observations to keep (1-alpha)
  if (pa$opt=="MIXT") tc.set=(rank(pre.z_m, ties.method="random")> floor(pa$n*(pa$alpha)))
  if (pa$opt=="HARD")    tc.set=(rank(pre.z_h, ties.method="random")> floor(pa$n*(pa$alpha)))
  
  ##### Assign points to best cluster, set trimmed points to cluster 0
  iter$assig <- apply(ll,1,which.max)*tc.set
  
  ##### POSTERIOR MEMBERSHIP (z_ij)
  if (pa$opt=="MIXT") iter$z_ij=ll/(pre.z_ + (pre.z_ ==0))*tc.set   #MIXT assigment including trimming 
  if (pa$opt=="HARD"){
    iter$z_ij=0*iter$z_ij
    ##### Set 1 for the assigned cluster, 0 otherwise
    iter$z_ij[cbind((1:pa$n),iter$assig+(iter$assig==0))] <- 1	
    iter$z_ij[tc.set==FALSE,] <- 0
  }
  
  #####	Convergence check for HARD assignment: check if assignments stabilized
  if (pa$opt=="HARD")      iter$code <- 2 * all(old.assig == iter$assig)		
  
  ##### UPDATE CLUSTER SIZES AND WEIGHTS
  if (pa$opt=="HARD")	  iter$csize <- tabulate (iter$assig, pa$K)
  if (pa$opt=="MIXT") iter$csize <- apply(iter$z_ij, 2,sum)  
  if (!pa$equal.weights)  iter$cw <- iter$csize/sum(iter$csize)	
  
  return (iter)
}

######## FUNCTION FOR obtaining the objective functions value for MIXT (obj_m) 
######## and HARD (obj_h) 

calcobj_tHDDC <- function (X, iter, pa) 
{			
  ww_m=matrix(0,nrow=pa$n,ncol=1)
  ww_h=matrix(0,nrow=pa$n,ncol=1)
  
  for (k in 1:pa$K) {
    A <- A.fun(X, pa, iter,k)
    w_m <- log(iter$cw[k])- 0.5*(A$A1+A$A2+A$A3+A$A4+A$A5)
    w_m[w_m==-Inf]=1e-99 
    ww_m=exp(w_m)+ww_m     
    w_h=w_m*(iter$assig==k)
    ww_h=w_h+ww_h     
  }
  if (pa$opt=="MIXT")		  iter$obj <- sum(log(ww_m[iter$assig>0]))  
  if (pa$opt=="HARD")			iter$obj <- sum(ww_h[iter$assig>0])  
  return (iter)
}


######## FUNCTION FOR obtaining the bic value
# Calculates BIC penalizing by the number of high-dimensional parameters

bic <- function(pa, iter){
  ##### Complexity includes: centers, weights, eigenvectors (Stiefel manifold), and eigenvalues (intrinsic + noise)
  ##### Also adjusted for constraints via restr.fact
  iter$bic_val=-2*iter$obj + ( pa$K*pa$p + pa$K - 1 + pa$p*sum(iter$q) - 
                                 0.5*sum(iter$q*(iter$q + 1)) + (sum(iter$q) - 1)*(1 - 1/pa$restr.fact[1]) + 1 
                               + (pa$K - 1)*(1 - 1/pa$restr.fact[2]) + 1 )*log(pa$n*(1-pa$alpha))
  return (iter)
}


######## AUXILIAR FUNCTION TO perfom the concentration steps max.iter times
# Performs the concentration steps until convergence or max.iter
tHDDC_update <- function(X, pa, iter, max.iter){
  
  for (i in 0:max.iter)
  {
    ##### Step 1: Enforce eigenvalue constraints on scatter structure
    iter <- f.restr_tHDDC(iter=iter,pa=pa)	
    
    #####  Step 2: Determine current assignment and apply trimming (MIXT models and HARD)
    iter <- findClustAssig_tHDDC(X,iter,pa)	 	
    
    ##### Check for early convergence (HARD assignment only)
    if ((iter$code==2) ||		##	if findClustAssig_tHDDC returned code=2 (convergence)
        (i == max.iter))		##	or we're in the last concentration step:
      break			##	break the for - loop - we finished this iteration!
    ##  don't re-estimate cluster parameters this time
    
    ##### Step 3: Re-estimate cluster parameters (M-step)
    iter <- estimClustPar_tHDDC(X, iter, pa)		## estimates the cluster's parameters 
  }
  
  iter <- calcobj_tHDDC (X, iter, pa)			## Final objective value
  iter <- bic (pa, iter)                  ## Final BIC
  
  return(iter)
}



################################################################################
############################# MAIN FUNCTION ####################################
################################################################################

# Required packages for parallel execution and reproducible RNG
# install.packages("foreach")
# install.packages("iterators")
library(doRNG)
library(doParallel)

set.seed(8128) #For reproducibility

## List of auxiliary functions to be exported to parallel workers
f_aux <- c( "catell", "restr2_eigenv_tHDDC", "adjust_eigenvalues","f.restr_tHDDC", "getini_tHDDC", "InitClusters_tHDDC", 
            "estimClustPar_tHDDC", "A.fun", "findClustAssig_tHDDC", "calcobj_tHDDC", "bic", 
            "tHDDC_update")

## tHDDC Parallelized Main Function
tHDDC_paralelizado_foreach <- function(X,q=NA,alpha = 0.05, nstart = 40, cstep1=2, cstep2=20, nkeep=5,  
                                       equal.weights = FALSE, restr.fact=5,zero.tol = 1e-16,  
                                       trace = 0,   opt="HARD", sol_ini_p = FALSE, sol_ini = NA, 
                                       catell_Y=TRUE,  catell_threshold=0.1, dim_max=20)  
{
  # --- Data Validation ---
  if (!is.numeric (X)) 
    stop ("parameter x: numeric matrix/vector expected")
  if( !is.matrix (X))
    X <- matrix (X, ncol = 1)
  
  n <- nrow (X)
  p <- ncol (X)
  no.trim <- floor(n*(1-alpha))
  K = length(q)  
  
  # Ensure nkeep does not exceed nstart
  nkeep <- min(nkeep, nstart)
  
  # --- Input Parameters List (pa) ---
  pa <- list (		
    n = n,			    
    p = p,	         
    alpha=alpha,		
    trimm = n-no.trim,	  
    no.trim = no.trim,	
    q = q,		
    K = K,		
    equal.weights = equal.weights,		
    zero.tol = zero.tol,		
    restr.fact=restr.fact,          
    opt = opt,              
    sol_ini_p = sol_ini_p, 
    ###################
    cstep1=cstep1,                  
    cstep2=cstep2,
    nstart=nstart,
    nkeep=nkeep,
    ###################
    catell_Y=catell_Y,
    catell_threshold=catell_threshold,
    dim_max=dim_max
  )
  
  
  # --- Internal State List (iter) ---
  iter <- list (			
    obj = NA,
    bic_val = NA,
    assig = array(0,n),			
    csize = array(NA,K),		 
    cw = rep(NA,K),		
    evalues = array(NA,c(p,K)),    
    evectors = list(),  
    center = array(NA,c(K,p)),	
    code = NA,		
    z_ij = matrix(0,nrow=n,ncol=K)	            
  )
  
  # --- Parallel Backend Initialization ---
  # Adjust threads for OpenMP if necessary and setup the cluster
  Sys.setenv(OMP_NUM_THREADS="4")
  parclus <- makeCluster(20)
  registerDoParallel(parclus)
  on.exit(stopCluster(parclus)) # Ensure cluster stops even if function fails
  
  
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
    iniciales <- foreach(inicial=1:pa$nstart, .export=f_aux, .inorder=TRUE) %dorng% {
      return(tHDDC_update(X, pa, InitClusters_tHDDC(X,iter,pa),pa$cstep1))}
  }
  
  
  print("j")  
  
  ##############################################################################
  ############# PHASE 2: SELECT THE nkeep BEST SOLUTIONS (BY BIC) ##############
  ##############################################################################
  
  ##### Sort solutions based on BIC (lower is better) and pick top candidates
  mejores <- order(sapply(iniciales, function(iter) iter$bic_val), decreasing=FALSE)[1:pa$nkeep]  
  
  print(mejores)
  
  ##############################################################################
  ############# PHASE 3: ITERATE TOP SOLUTIONS UNTIL CONVERGENCE ###############
  ##############################################################################
  
  ##### Run the selected best candidates for a longer duration (cstep2)
  finales <- foreach(l=mejores, .export=f_aux) %dorng% {
    return(tHDDC_update(X, pa, iniciales[[l]], pa$cstep2))} 
  
  ##### Find the overall best solution from the refined results
  j <- which.min(sapply(finales, function(iter) iter$bic_val))
  
  ##### Final output
  return(finales[[j]])
}  

