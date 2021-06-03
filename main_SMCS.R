# Script allowing to perform robust Bayesian optimization with scarlar modeling of constraints and selection of most useful constraint

library('DiceDesign')
library('DiceKriging')
library("parallel")
library("Rsolnp")
library("nloptr")
library("randtoolbox")
library('reshape2')
library('future.apply')
library('doFuture')

########################################
######################################## krigeage sur l ensemble des donnees 
######################################## + creation des fonctions objective/contraintes a partir des 900 points
######################################## + R_functions


mean_sd_Z <- function(x,model,alea,d,m,select){
  x <- matrix(x,ncol=d,nrow=1)
  alea <- matrix(alea,ncol=m)
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  pred <- DiceKriging::predict.km(model,dat,checkNames = FALSE,type="SK",cov.compute=FALSE)
  if(select=="TRUE") return(mean(pred$mean))
  if(select=="FALSE") return(abs(mean(pred$sd)))
}

Phi_G <- function(x,modelConst,alea,seuil,d,m,sign){
  input <- matrix(0,dim(alea)[1],d+m)
  input[,((d+1):(d+m))] <- alea
  for(j in 1:dim(alea)[1]) input[j,1:d] <- x
  pred <- DiceKriging::predict.km(modelConst,input,checkNames = FALSE,type="SK",cov.compute=FALSE)
  if(sign == ">") phi <- pnorm((pred$mean-seuil)/pred$sd)
  else phi <- pnorm((seuil-pred$mean)/pred$sd)
  return(phi)
}

Expectation_A <- function(x,alea,d,m,list_const_km,seuil,sign,alpha){
  d_c <- length(list_const_km) # constraints models
  PHI <- matrix(NA,nrow=d_c,ncol=dim(alea)[1])
  for(i in 1:d_c) PHI[i,] <- Phi_G(x,list_const_km[[i]],alea,seuil[i],d,m,sign[i])
  E_A <- 1 - alpha - mean(apply(PHI,2,prod))
  return(E_A)
}

z_mean_feasible <- function(modelObjective,d,m,alea,list_const_km,seuil,sign,alpha){
  X <- modelObjective@X[,1:d]
  res <- apply(X,1,mean_sd_Z,modelObjective,alea,d,m,select="TRUE")
  res1 <- apply(X,1,Expectation_A,alea,d,m,list_const_km,seuil,sign,alpha)
  index <- which(res1 <= 0)
  if(length(index) > 0 ) {index2 <- which.min(res[index]);   return(list(min=res[index[index2]],opt=index[index2],mar=1))}
  if(length(index) == 0 ) {t <- which.min(res1); return(list(min=res[t],opt=t,mar=2))}
}

Probability_A <- function(x,alea,Nsim,d,m,list_const_km,seuil,sign,alpha){
  input <- matrix(0,dim(alea)[1],d+m)
  input[,(d+1):(d+m)] <- alea
  for(j in 1:dim(alea)[1]) input[j,1:d] <- x
  d_c <- length(list_const_km) # constraints models
  realizations <- array(0, dim=c(Nsim,dim(alea)[1],d_c) ) # multidimensional array (Nsim,CRN,nbreConstraint)
  for(i in 1:d_c) realizations[,,i] <- DiceKriging::simulate(list_const_km[[i]],nsim=Nsim,newdata=input,
                                                             checkNames = FALSE,type="SK",cond=TRUE)
  for(i in 1:d_c)
  {
    if(sign[i] == ">") realizations[,,i] <-  seuil[i] - realizations[,,i]
    else realizations[,,i] <- realizations[,,i] - seuil[i]
  }
  sim <- apply(realizations,c(1,2),max) ## max(g1,g2,g3)
  p_sim <- apply(sim,1,function(vect) length(which(vect<=0))/dim(alea)[1])
  c <- length(which(1 - alpha - p_sim <= 0))/Nsim
  return(c=c)
}

Feas_Expected_Improvement <- function(x,minimum_feas,model,alea,Nsim,d,m,list_const_km,seuil,sign,alpha){
  sdZ <- abs(mean_sd_Z(x,model,alea,d,m,select="FALSE"))
  espZ <- mean_sd_Z(x,model,alea,d,m,select="TRUE")
  v <- (minimum_feas-espZ)/sdZ
  phi <- dnorm(v, mean = 0, sd = 1, log = FALSE)
  PHI <- pnorm(v, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)  
  EIZ <- sdZ*(v*PHI+phi)
  proba <- Probability_A(x,alea,Nsim,d,m,list_const_km,seuil,sign,alpha)
  return(c(EIZ)*proba)
}

PredVar_Z <- function(xstar,xp1,model,alea,d,m){
  U_MC <- dim(alea)[1]
  xxstar <- xstar
  objectinit <- model@covariance
  
  X1 <- (as.matrix( rbind(model@X,as.vector(xp1))))
  X2 <- (as.matrix( rbind(model@X,as.vector(xp1))))
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  dim <- ncol(X1)
  # C <- covMat1Mat2(model@covariance,X1,X2,nugget.flag = TRUE)/model@covariance@sd2 
  C <- covMatrix(model@covariance,X1, noise.var = 1e-8)$C/model@covariance@sd2
  C <- solve(C)
  
  ########################################################### FIRST TERM
  U1 <- alea
  U2 <- alea
  nU1 <- nrow(U1)
  nU2 <- nrow(U2)
  dimU <- ncol(U1)
  object <- objectinit
  object@range.val <- covparam2vect(objectinit)[(d+1):(d+m)]
  outU <- covMat1Mat2(object,U1,U2,nugget.flag = TRUE)/model@covariance@sd2 
  MU <- mean(matrix(outU, nU1, nU2))
  term1 <- MU
  
  ########################################################### SECOND TERM
  
  X1 <- cbind(t(matrix(rep(xstar,U_MC),length(xstar),U_MC)),alea)
  X2 <- (as.matrix( rbind(model@X,as.vector(xp1))))
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  dim <- ncol(X1)
  object <- objectinit
  object@range.val <- covparam2vect(objectinit)
  M <- covMat1Mat2(model@covariance,X1,X2,nugget.flag = TRUE)/model@covariance@sd2
  
  X1 <- (as.matrix( rbind(model@X,as.vector(xp1))))
  X2 <- cbind(t(matrix(rep(xstar,U_MC),length(xstar),U_MC)),alea) 
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  dim <- ncol(X1)
  object <- objectinit
  object@range.val <- covparam2vect(objectinit)
  trans_M <- covMat1Mat2(model@covariance,X1,X2,nugget.flag = TRUE)/model@covariance@sd2
  
  a <- matrix(apply(M,2,mean),nrow=1)
  b <- matrix(apply(trans_M,1,mean),ncol=1)
  term2 <- a%*%C%*%(b)
  ###########################################################   
  return(term1 - term2)
}  

PredMean_Z <- function(x,xp1,observation,model,alea,d,m){
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  colnames(dat) <- NULL
  xp1 <- as.vector(xp1)
  data <- data.frame(rbind(xp1,dat))
  utile <- DiceKriging::predict.km(model,data,checkNames = FALSE,type="SK",cov.compute=TRUE)
  term1 <- mean_sd_Z(x,model,alea,d,m,"TRUE")
  term2 <- (observation - utile$mean[1])/utile$sd[1]^2
  term3 <- mean(utile$cov[1,])  
  return(term1 + term2*term3)
} 

PredVar_F <- function(newdata,xp1,model,d,m){
  newdata <- matrix(newdata,ncol=d+m)
  xp1 <- matrix(xp1,ncol=d+m)
  utile <- DiceKriging::predict.km(model,rbind(xp1,newdata),checkNames = FALSE,type="SK",cov.compute=TRUE)
  return(utile$sd[2]^2 - utile$cov[1,2]^2/utile$sd[1]^2)
}

TotalVar <- function(x,xp1,model,alea,d,m,observation,ming,variZplus1){
  meanZplus1 <- PredMean_Z(x=x,xp1=xp1,observation,model=model,alea=alea,d=d,m=m)
  
  term <- (ming - meanZplus1)/sqrt(variZplus1)
  a <- pnorm(term, mean = 0, sd = 1)
  b <- dnorm(term, mean = 0, sd = 1)
  
  term1 <- (ming - meanZplus1)^2 + variZplus1
  term2 <- sqrt(variZplus1)*(ming - meanZplus1)
  
  expected_imprmnt <- (ming - meanZplus1)*a + sqrt(variZplus1)*b
  
  return(c(expected_imprmnt,term1*a + term2*b - expected_imprmnt^2))
}


u_tp1_f <- function(utilde,xnew,rep,ming,modelObjective,alea,m,d){
  input4s <- c(xnew,utilde)
  stopping <- 0
  oldlaw <- DiceKriging::predict.km(modelObjective,rbind(input4s),checkNames = FALSE,type="SK",cov.compute=TRUE)
  oldlaw <- c(oldlaw$mean,oldlaw$sd)
  aleaa <- matrix(randtoolbox::sobol(n=rep, dim = 1, init = TRUE, scrambling = 0, seed = 4711, normal = TRUE),ncol=1)
  aleaa <- aleaa*oldlaw[2] + oldlaw[1]
  variZplus1 <- PredVar_Z(xnew,input4s,modelObjective,alea,d,m)
  term1 <- term2 <- NULL
  repeat{
    stopping = stopping+1
    muet <- TotalVar(xnew,input4s,modelObjective,alea,d,m,aleaa[stopping],ming,variZplus1)
    term1 <- c(term1,muet[1]) ;   term2 <- c(term2,muet[2])
    if (stopping == rep) break
  }
  res <- var(term1) + mean(term2)
  return(res)
}


u_tp1_g <- function(utilde,xnew,inputsss,m,d,list_const_km,seuil,sign, const_ind){
  input4s <- c(xnew,utilde)
  d_c <- length(list_const_km) # constraints models
  PHI <- matrix(NA,nrow=d_c,ncol=dim(inputsss)[1])
  for(i in 1:d_c)
  {
    pred <- DiceKriging::predict(list_const_km[[i]],newdata=inputsss,checkNames = FALSE,type="SK",cond=TRUE)
    if(i == const_ind) predicted_sd2 <- apply(inputsss,1,PredVar_F,input4s,list_const_km[[i]],d,m)
    else predicted_sd2 <- pred$sd^2
    if(sign[i] == ">") PHI[i,] <- pnorm((pred$mean-seuil[i])/sqrt(abs(predicted_sd2)))
    else PHI[i,] <- pnorm((seuil[i]-pred$mean)/sqrt(abs(predicted_sd2)))
  }
  pq <- apply(PHI,2,prod)
  pp <- mean(pq*(1-pq))
  return(pp)
}

########################################
######################################## INITIAL DESIGN OF EXPERIMENTS & metamodels
######################################## 


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


d <- 2 # Dimension of X
m <- 2 # Dimension of U        
n_g = 2 # Number of constraints
n = 6*(d+m) #Number of points in the inital data set
U_MC <- 300 # Montecarlo generation of samples of U
N_xOpt <- 300 # Search grid for the optimal values of X
N_uOpt <- 300 # Search grid for the optimal values of U


Nsim <- 100 # Number of simulated trajectories for the computation of P(C < 0) 
alpha <- 0.05 # Tollerated probability of failure
beta <- 1 - alpha
quantizer <- 10

seuil <- c(0.,0.) # threshold
sign <- c("<","<")
InpNames = c('x1','x2','u1','u2')
ConstNames = c('g1','g2')


iterations = 80 #Infilled data samples for the optimization


################ 
################ Initial Doe and simulator responses
################ 

# Sampling of the initial data set. For this test-case, two independent uniform distributions are considered for u1 and u2.
design = matrix(maximinSA_LHS(lhsDesign(n,d, seed = 0)$design)$design, ncol = d)
design = cbind(design,runif(n))
design = cbind(design,runif(n))
design_init = design 

alea <- matrix(cbind(runif(U_MC),runif(U_MC)),ncol = m)

lhs_xtarg = lhsDesign(N_xOpt,d)$design

Opt = lhsDesign(N_uOpt,m)$design

#simulator responses
response_f <- apply(design, 1, Objective)
response_g_1 <- apply(design, 1, Constraint1)
response_g_2 <- apply(design, 1, Constraint2)

df = data.frame(design,response_g_1,response_g_2)
df = setNames(df,c(InpNames,ConstNames))

df = melt(df, id.vars = InpNames,variable.name = "output",
          value.name = "value")

colnames(design) = InpNames
  
  

response_f <- apply(design_init, 1, Objective) #simulator responses
response_g_1 <- apply(design_init, 1, Constraint1) #simulator responses
response_g_2 <- apply(design_init, 1, Constraint2) #simulator responses

covtype <- "matern5_2"

model_O <- km(formula= ~1, design = data.frame(x=design_init), response= response_f,covtype=covtype, nugget=1e-5) # Objective function metamodel (F)
model_C1 <- km(formula= ~1, design = data.frame(x=design_init), response= response_g_1,covtype=covtype, nugget=1e-5) # Constraint function metamodel (G)
model_C2 <- km(formula= ~1, design = data.frame(x=design_init), response= response_g_2,covtype=covtype, nugget=1e-5) # Constraint function metamodel (G)


listCon <- list(model_C1=model_C1,model_C2=model_C2)

min_feas_init <- z_mean_feasible(model_O,d,m,alea,listCon,seuil,sign,alpha) # feasible minimum (z_min^{feas})
hist_feasmin <- c(min_feas_init$min)  # feasible minimum value 
hist_feasopt <- c(min_feas_init$opt) # associated point 

print("initialization: Done")

########################################
######################################## Optimization iterations
######################################## 
design_init_g1 = design_init
design_init_g2 = design_init

for(i in 1:iterations)
{
  
  cat('Iteration number ', i)
  
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  EFI <- future_apply(lhs_xtarg,1,Feas_Expected_Improvement,min_feas_init$min,model_O,alea,Nsim,d,m,listCon,seuil,sign,alpha)
  x_tplus1 <- lhs_xtarg[which.max(EFI),]
  
  ##### bobyqa in U
  
  Opt_ <- matrix(0,N_uOpt,(d+m))
  Opt_[,(d+1):(d+m)] <- Opt
  for(j in 1:N_uOpt) Opt_[j,1:d] <- x_tplus1

  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  optuf = future_apply(Opt,1,u_tp1_f,x_tplus1,rep=5,min_feas_init$min,model_O,alea,m,d)
  u_tplus1f <- Opt_[which.min(optuf),(d+1):(d+m)]
  
  u_tplus1g_list = list()
  varval_list = c()
  for(const_ind in 1:length(listCon))
  {
    registerDoFuture()
    plan(multisession, workers = max(detectCores()-1, 1))
    optug = future_apply(Opt,1,u_tp1_g,x_tplus1,Opt_,m,d,listCon,seuil,sign, const_ind)
    optug[which(is.nan(optug) == TRUE)] = 1e3
    u_tplus1g_list[[const_ind]] = Opt_[which.min(optug),(d+1):(d+m)]
    varval_list = c(varval_list, min(optug))
  }
  
  select_const = which.min(varval_list)
  u_tplus1g <- u_tplus1g_list[[select_const]]

      ##### Update metamodel
  
  point_new1 = cbind(matrix(x_tplus1,1,d),matrix(u_tplus1f,1,m))
  point_new2 = cbind(matrix(x_tplus1,1,d),matrix(u_tplus1g,1,m))
  
  design_init = rbind(design_init, point_new1)
  response_f = c(response_f, Objective(point_new1))
  model_O <- km(formula= ~1, design = data.frame(x=design_init), response= response_f,covtype=covtype, nugget=1e-5) # Objective function metamodel (F)
  
  
  if(select_const == 1)
  {
    response_g_1 = c(response_g_1, Constraint1(point_new2)) 
    design_init_g1 = rbind(design_init_g1, point_new2)
    model_C1 <- km(formula= ~1, design = data.frame(x=design_init_g1), response= response_g_1,covtype=covtype, nugget=1e-5) # Constraint function metamodel (G)
  }

  if(select_const == 2)
  {
    response_g_2 = c(response_g_2, Constraint2(point_new2)) 
    design_init_g2 = rbind(design_init_g2, point_new2)
    model_C2 <- km(formula= ~1, design = data.frame(x=design_init_g2), response= response_g_2,covtype=covtype, nugget=1e-5) # Constraint function metamodel (G)
  }
  
  print(paste("metamodel_update ",i," : Done",sep=""))
  
  listCon <- list(model_C1=model_C1,model_C2=model_C2)
  
  min_feas_init <- z_mean_feasible(model_O,d,m,alea,listCon,seuil,sign,alpha) # feasible minimum (z_min^{feas})
  hist_feasmin <- c(hist_feasmin,min_feas_init$min)  # feasible minimum value 
  hist_feasopt <- c(hist_feasopt,min_feas_init$opt) # associated point 
  
}


design=design_init

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



