source("Auxiliary_functions.R")

# Definitions of the problem functions
Constraint1 <- function(X)
{
  x1 = X[1]*10-5
  x2 = X[2]*10-5
  u1 = X[3]*10-5
  u2 = X[4]*10-5
  
  g1 = -x1^2 + 5*x2 - u1 + u2^2 -1
  
  return(g1/50)
}

Constraint2 <- function(X)
{
  x1 = X[1]*10-5
  x2 = X[2]*10-5
  u1 = X[3]*10-5
  u2 = X[4]*10-5
  
  g2 = (-x1^2 + 5*x2 - u1 + u2^2 -1)*(x1+5)/5 - u1 - 1
  
  return(g2/30)
}

Objective <- function(X)
{
  x1 = X[1]*10-5
  x2 = X[2]*10-5
  u1 = X[3]*10-5
  u2 = X[4]*10-5
  
  f = 5*(x1^2+x2^2) - (u1^2+u2^2) + x1*(u2-u1+5) + x2*(u1-u2+3)
  return(f)
}

# Definitions of a single function allowing to call either constraints depending on the value of an additional variable x5
constraintsFun = function(x)
{
  if(x[5] == 'g1') return( Constraint1(x[1:(d+m)]))
  
  else if (x[5] == 'g2') return(Constraint2(x[1:(d+m)]))
  
  else print('Error')
}


d <- 2 # Dimension of X
m <- 2 # Dimension of U        
n_g = 2 # Number of constraints
n = 6*(d+m) #Number of points in the inital data set
U_MC <- 300 # Montecarlo generation of samples of U
N_xOpt <- 2000 # Search grid for the optimal values of X
N_uOpt <- 1000 # Search grid for the optimal values of U


Nsim <- 100 # Number of simulated trajectories for the computation of P(C < 0) 
seuil <- 0. # Limit value of the constraints
alpha <- 0.05 # Tollerated probability of failure
beta <- 1 - alpha
quantizer <- 10


iterations = 160 #Infilled data samples for the optimization
InpNames = c('x1','x2','u1','u2')
ConstNames = c('g1','g2')

################ 
################ Initial Doe and simulator responses
################ 

# Sampling of the initial data set. For this test-case, two independent uniform distributions are considered for u1 and u2.
design = matrix(maximinSA_LHS(lhsDesign(n,d, seed = rep)$design)$design, ncol = d)
design = cbind(design,runif(n))
design = cbind(design,runif(n))

alea <- matrix(cbind(runif(U_MC),runif(U_MC)),ncol = m)

lhs_xtarg = lhsDesign(N_xOpt,d)$design

Opt = lhsDesign(N_uOpt,m)$design

#simulator responses
responseObjective <- apply(design, 1, Objective)
responseConstraint1 <- apply(design, 1, Constraint1)
responseConstraint2 <- apply(design, 1, Constraint2)

colnames(design) = InpNames

#Type of covariances used for the continuous variables
covtype <- "matern5_2"

for(iter in 1:iterations){
  ################ 
  ################ Initial GP models
  ################ 
  repeat{
    modelObjective <- km(formula= ~1, design = data.frame(x=design), response=responseObjective, 
                         covtype=covtype, nugget=1e-5)
    
    if(min(round(modelObjective@covariance@range.val,2)) != 0 &&  modelObjective@covariance@sd2 != 0) break
  }
  
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  multistart <- 35
  
  #Continuous kernels 
  kContx1 <- covRadial(k1Fun1 = k1Fun1Matern5_2,
                       d = 1, cov = "corr", input = "x1")
  coefUpper(kContx1) <- 100
  coefLower(kContx1) <- 0.
  
  kContx2 <- covRadial(k1Fun1 = k1Fun1Matern5_2,
                       d = 1, cov = "corr", input = "x2")
  coefUpper(kContx2) <- 100
  coefLower(kContx2) <- 0.
  
  kContu1 <- covRadial(k1Fun1 = k1Fun1Matern5_2,
                       d = 1, cov = "corr", input = "u1")
  coefUpper(kContu1) <- 100
  coefLower(kContu1) <- 0.
  
  kContu2 <- covRadial(k1Fun1 = k1Fun1Matern5_2,
                       d = 1, cov = "corr", input = "u2")
  coefUpper(kContu2) <- 100
  coefLower(kContu2) <- 0.
  
  # Discrete kernel
  kCatout <- q1Symm(factor = df$output, input = "output", cov = "homo")
  coefUpper(kCatout) <- c(pi,300.)
  coefLower(kCatout) <- c(0., 1e-3)
  
  # Composed kernel
  kMix <- covComp(formula = ~ kContx1() * kContx2() * kContu1() * kContu2() * kCatout() )
  
  
  opts <- list(ftol_abs = 1e-7, maxeval = 10000)
  
  # Multiple-output surrogate model of the constraints
  modelConstraint <-  gp(formula = value ~ 1 ,
                         cov = kMix, data = df, multistart = multistart,  opts = opts ,
                         optimMethod = "NLOPT_LN_COBYLA",  compGrad = FALSE, trace = 0, varNoiseLower = 1e-6, varNoiseUpper = 1e-3) 
  
  
  ################
  ################ Current Feasible Minimum
  
  #############  ################
  hist_feasminf <- hist_feasopt <- NULL
  feasminf <- minimin_g(modelObjective,modelConstraint,d,m,alea,beta,n_g)
  hist_feasminf <- c(hist_feasminf,feasminf$min)
  hist_feasopt <- rbind(feasminf$opt) # associated point
  ###
  ################ Find x_{t+1} by maximizing EFI ("discrete" solver)
  ################
  
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  EFI <- future_apply(lhs_xtarg,1,Feas_Expected_Improvement,feasminf$min,modelObjective,alea,d,m,Nsim,modelConstraint,beta,n_g)
  x_tplus1 <- lhs_xtarg[which.max(EFI),]
  
  #   ################ 
  #   ################ Find u_{t+1} by minimizing S ("Discrete Solver)
  #   ################ 
  
  Opt_ <- matrix(0,N_uOpt,(d+m))
  Opt_[,(d+1):(d+m)] <- Opt
  for(j in 1:N_uOpt) Opt_[j,1:d] <- x_tplus1
  newxu<- Sampling_new_singleg(x_tplus1,rep=quantizer,modelObjective,alea,feasminf$min,Opt_,modelConstraint,m,d,n_g)
  
  # Data set is updated
  newxu_f = newxu[[1]]
  newxu_g1 = newxu[[2]] 
  
  newvalue1 = constraintsFun(newxu_g1)
  
  newsample1 = cbind(newxu_g1,newvalue1)
  newsample1 = setNames(newsample1,names(df))
  
  
  df = rbind(df,newsample1)
  design = rbind(design,newxu_f[1:(d+m)])
  newObj  = Objective(newxu_f[1:(d+m)])
  responseObjective = c(responseObjective,newObj[1])
  
  res = list(df,design)
  RES[[rep]] = res
  
  fichier_res <- "Res_MMCS"
  ff <- paste(fichier_res, "RData", sep=".")
  save(RES,file=ff)
  }







