# Script allowing to perform robust Bayesian optimization with multi-output modeling of constraints

source("Auxiliary_functions.R")

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
N_xOpt <- 300 # Search grid for the optimal values of X
N_uOpt <- 300 # Search grid for the optimal values of U


Nsim <- 100 # Number of simulated trajectories for the computation of P(C < 0) 
seuil <- 0. # Limit value of the constraints
alpha <- 0.05 # Tollerated probability of failure
beta <- 1 - alpha
quantizer <- 10


iterations = 10 #Infilled data samples for the optimization
InpNames = c('x1','x2','u1','u2')
ConstNames = c('g1','g2')

################ 
################ Initial Doe and simulator responses
################ 

# Sampling of the initial data set. For this test-case, two independent uniform distributions are considered for u1 and u2.
design = matrix(maximinSA_LHS(lhsDesign(n,d, seed = 0)$design)$design, ncol = d)
design = cbind(design,runif(n))
design = cbind(design,runif(n))

alea <- matrix(cbind(runif(U_MC),runif(U_MC)),ncol = m)

lhs_xtarg = lhsDesign(N_xOpt,d)$design

Opt = lhsDesign(N_uOpt,m)$design

#simulator responses
responseObjective <- apply(design, 1, Objective)
responseConstraint1 <- apply(design, 1, Constraint1)
responseConstraint2 <- apply(design, 1, Constraint2)

df = data.frame(design,responseConstraint1,responseConstraint2)
df = setNames(df,c(InpNames,ConstNames))

df = melt(df, id.vars = InpNames,variable.name = "output",
          value.name = "value")

colnames(design) = InpNames
  
covtype <- "matern5_2"

for(iter in 1:iterations){
  
  cat('Iteration number ', iter)

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
  multistart <- 30
  
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
  
  kCatout <- q1Symm(factor = df$output, input = "output", cov = "homo")
  coefUpper(kCatout) <- c(pi,300.)
  coefLower(kCatout) <- c(0., 1e-3)
  
  kMix <- covComp(formula = ~ kContx1() * kContx2() * kContu1() * kContu2() * kCatout() )
  
  
  opts <- list(ftol_abs = 1e-7, maxeval = 10000)
  
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
  ################ Find x_{t+1} by EFI ("discrete" or "bobyqa" solver)
  ################
  
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  EFI <- future_apply(lhs_xtarg,1,Feas_Expected_Improvement,feasminf$min,modelObjective,alea,d,m,Nsim,modelConstraint,beta,n_g)
  x_tplus1 <- lhs_xtarg[which.max(EFI),]
  
  #   ################ 
  #   ################ Find u_{t+1} by S ("Discrete" or "Bobyqa" Solver)
  #   ################ 
  
  Opt_ <- matrix(0,N_uOpt,(d+m))
  Opt_[,(d+1):(d+m)] <- Opt
  for(j in 1:N_uOpt) Opt_[j,1:d] <- x_tplus1
  newxu<- Sampling_new_singleu(x_tplus1,rep=quantizer,modelObjective,alea,feasminf$min,Opt_,modelConstraint,m,d,n_g)
  
  newxu_f = newxu[[1]]
  newxu_g1 = newxu[[1]]
  newxu_g2 = newxu[[2]]
  
  newvalue1 = Constraint1(newxu_g1)
  newvalue2 = Constraint2(newxu_g2)
  
  newsample1 = cbind(newxu_g1,newvalue1)
  newsample1 = setNames(newsample1,names(df))
  
  newsample2 = cbind(newxu_g2,newvalue2)
  newsample2 = setNames(newsample2,names(df))
  
  df = rbind(df,newsample1,newsample2)
  design = rbind(design,newxu_f[1:(d+m)])
  newObj  = Objective(newxu_f[1:(d+m)])
  responseObjective = c(responseObjective,newObj[1,1])
  
}



#Post Processing results

faux = function(X)
{
  
  Xl = matrix(0, nrow = 0, ncol = d)
  for(j in 1:1000) Xl = rbind(Xl,X)
  Inp = cbind(Xl,U)
  
  F = apply(Inp,1,Objective)
  
  G1 = apply(Inp,1,Constraint1)
  
  G2 = apply(Inp,1,Constraint2)
  
  
  return(c(mean(F),length(which(G1<= 0. & G2 <= 0))/1000))
}

#Montecarlo runs used to computed the expected value of the objective function and the probability of feasibility on the sampled data points
U = matrix(cbind(runif(1000),runif(1000)),ncol = m) 
registerDoFuture()
plan(multisession, workers = max(detectCores()-1, 1))

X =design[,1:(d)]

registerDoFuture()
plan(multisession, workers = max(detectCores()-1, 1))
Processed_res = future_apply(X,1,faux)


fmean = Processed_res[1,]
pg =  Processed_res[2,]

Optimum_value = min(fmean[which(pg>=beta)])

cat('Identified optimal robust value equal to ', Optimum_value)



