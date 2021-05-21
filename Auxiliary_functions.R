library('DiceDesign')
library('DiceKriging')
library("parallel")
library("Rsolnp")
library("nloptr")
library('mvtnorm')
library('kergp')
library('reshape2')
library('doFuture')
library('future.apply')
source("simulate_new.R")


###################################################################################### 
########################################### Function to perform quantization
###################################################################################### 
choice.grid <- function(X, N, ng = 1, p = 2) {
  if (!is.numeric(X)) 
    stop("X must be numeric")
  if (!is.vector(X) & !is.matrix(X)) 
    stop("X must be a matrix or a vector")
  if ((!(floor(N) == N)) | (N <= 0)) 
    stop("N must be entire and positive")
  if ((!(floor(ng) == ng)) | (ng <= 0)) 
    stop("B must be entire")
  if (p < 1) 
    stop("p must be at least 1")
  if (is.vector(X)) {
    d <- 1
    n <- length(X)
    primeX <- matrix(sample(X, n * ng, replace = TRUE), 
                     nrow = ng)
    hatX <- replicate(ng, sample(unique(X), N, replace = FALSE))
    hatX0 <- hatX  #save the initial grids
    gammaX <- array(0, dim = c(ng, n + 1))
    # initialisation of the step parameter gamma
    a_gamma <- 4 * N
    b_gamma <- pi^2 * N^(-2)
    BtestX = array(Inf, dim = c(ng, 1))
    # choice of gamma0X
    projXbootinit <- array(0, dim = c(n, ng))
    # index of the grid on which X is projected
    iminx <- array(0, dim = c(n, ng))
    for (i in 1:n) {
      RepX <- matrix(rep(X[i], N * ng), ncol = ng, byrow = TRUE)
      Ax <- sqrt((RepX - hatX0)^2)
      iminx[i, ] <- apply(Ax, 2, which.min)
      mx <- matrix(c(iminx[i, ], c(1:ng)), nrow = ng)
      projXbootinit[i, ] <- hatX0[mx]
    }
    RepX <- matrix(rep(X, ng), ncol = ng)
    distortion <- apply((RepX - projXbootinit)^2, 2, sum)/n
    temp_gammaX <- which(distortion > 1)
    if (length(temp_gammaX) > 0) {
      distortion[temp_gammaX] <- array(1, dim = c(length(temp_gammaX), 
                                                  1))
    }
    gamma0X <- distortion
    if(any(gamma0X < 0.005)){gamma0X[gamma0X<0.005] <- 1}
    gammaX = array(0, dim = c(ng, n + 1))
    for (i in 1:ng) {
      # calculation of the step parameter
      gammaX[i, ] <- gamma0X[i] * a_gamma/(a_gamma + gamma0X[i] * 
                                             b_gamma * c(1:(n + 1) - 1))
    }
    gammaX[, 1] <- gamma0X
    iminX <- array(0, dim = c(n, ng))  #index that will change 
    tildeX <- array(0, dim = c(N, ng))
    # updating of the grids, providing optimal grids
    for (i in 1:n) {
      for (j in 1:ng) {
        tildeX[, j] <- matrix(rep(primeX[j, i], N), nrow = N, 
                              byrow = TRUE)
      }
      Ax <- sqrt((tildeX - hatX)^2)
      # calculation of each distance to determine the point of
      # the grid the closer of the stimuli
      iminX[i, ] <- apply(Ax, 2, which.min)
      mX <- matrix(c(iminX[i, ], c(1:ng)), nrow = ng)
      if(sqrt(sum((hatX[mX] - primeX[, i])^2))==0){
        hatX[mX] <- hatX[mX] - gammaX[, i + 1] * (hatX[mX] - primeX[, i])
      }else{
        hatX[mX] <- hatX[mX] - gammaX[, i + 1] * (hatX[mX] - 
                                                    primeX[, i]) * (sqrt(sum((hatX[mX] - primeX[, 
                                                                                                i])^2)))^(p - 1)/sqrt(sum((hatX[mX] - primeX[,                                                                                                                          i])^2))
      }
    }
  } else {
    n <- ncol(X)
    d <- nrow(X)
    primeX <- array(X[, sample(c(1:n), n * ng, 
                               replace = T)], dim = c(d, n, ng))
    hatX <- replicate(ng, unique(X)[, sample(c(1:n), 
                                             N, replace = F)])  # initial grids chosen randomly in the sample
    hatX <- array(hatX, dim = c(d, N, ng))
    hatX0 <- hatX  #save the initial grids
    # initialisation of the step parameter gamma
    a_gamma <- 4 * N^(1/d)
    b_gamma <- pi^2 * N^(-2/d)
    BtestX <- array(Inf, dim = c(1, ng))
    # choice of gamma0X
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        Bx <- array(0, dim = c(ng, 1))
        Bx <- sqrt(apply((hatX[, i, , drop = FALSE] - 
                            hatX[, j, , drop = FALSE])^2, c(2, 3), sum))/2
        temp_gammaX <- which(Bx < BtestX)
        BtestX[temp_gammaX] <- Bx[temp_gammaX]
      }
    }
    temp_gammaX = which(BtestX > 1)
    if (length(temp_gammaX) > 0) {
      BtestX[temp_gammaX] <- array(1, dim = c(length(temp_gammaX), 
                                              1))
    }
    gamma0X <- BtestX
    if(any(gamma0X < 0.005)){gamma0X[gamma0X<0.005] <- 1}
    
    gammaX <- array(0, dim = c(ng, n + 1))
    for (i in 1:ng) {
      # calculation of the step parameter
      gammaX[i, ] <- gamma0X[i] * a_gamma/(a_gamma + gamma0X[i] * 
                                             b_gamma * c(1:(n + 1) - 1))
    }
    gammaX[, 1] <- gamma0X
    iminX <- array(0, dim = c(n, ng))  #index that will change 
    tildeX <- array(0, dim = c(d, N, ng))
    # updating of the grids, providing optimal grids
    for (i in 1:n) {
      for (j in 1:ng) {
        tildeX[, , j] <- matrix(rep(primeX[, i, j], N), 
                                nrow = N, byrow = FALSE)
      }
      Ax <- sqrt(apply((tildeX - hatX)^2, c(2, 3), sum))
      # calculation of each distance to determine the point of the grid 
      #the closer of the stimuli
      iminX[i, ] <- apply(Ax, 2, which.min)
      for (k in 1:d) {
        m <- matrix(c(rep(k, ng), iminX[i, ], c(1:ng)), ncol = 3)
        if(p==2){
          hatX[m] = hatX[m] - gammaX[, i + 1] * (hatX[m] - primeX[k, i, ]) 
        }else{
          hatX[m] = hatX[m] - gammaX[,i + 1] * (hatX[m] - primeX[k, i, ]) * 
            (sqrt(sum((hatX[m] - primeX[k, i, ])^2)))^(p - 1)/sqrt(sum((hatX[m] - primeX[k, 
                                                                                         i, ])^2))
        } 
      }
    }
  }
  output <- list(init_grid=hatX0,opti_grid=hatX)
  output
} 


###################################################################################### 
########################################### Integrated Process : mean & standard deviation
###################################################################################### 
IntegratedProcess <- function(x,model,alea,d,m,select)
{
  x <- matrix(x,ncol=d,nrow=1)
  alea <- matrix(alea,ncol=m)
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  pred <- DiceKriging::predict.km(model,dat,checkNames = FALSE,type="SK",cov.compute=FALSE)
  if(select=="TRUE") return(mean(pred$mean))
  if(select=="FALSE") return(abs(mean(pred$sd)))
}

###################################################################################### 
########################################### Expectation of the proba of violation
###################################################################################### 
Expectation_C <- function(x,modelConst,alea,d,m,n_g)
{
  
  input <- matrix(0,dim(alea)[1],d+m)
  input[,((d+1):(d+m))] <- alea
  for(j in 1:dim(alea)[1]) input[j,1:d] <- x 
  
  inp = data.frame()
  for (lv in levels(modelConst$X$output))
  {
    inp = rbind(inp,data.frame(input,factor(rep(lv,dim(alea)[1]))))
  }
  inp = setNames(inp, modelConst$inputNames)
  
  pred = predict(modelConst,inp,type="SK",cond=TRUE,covCompute = TRUE)
  
  toto = c()
  for(i in 1:dim(alea)[1])
  {
    m_l = c()
    for(j in 1:n_g)
    {
      m_l = c(m_l, pred$mean[(i+(j-1)*dim(alea)[1])])
    }
    
    cov_l = matrix(0, nrow = n_g, ncol = n_g)
    for(j in 1:n_g)
    {
      for(k in 1:n_g)
      {
        cov_l[j,k] = pred$cov[(i+(j-1)*dim(alea)[1]),(i+(k-1)*dim(alea)[1])]  
      }
    }
    mean = m_l
    sigma = cov_l
    lower <- rep(-Inf, n_g)
    upper <- rep(0, n_g)
    # lower <- rep(0, n_g)
    # upper <- rep(Inf, n_g)
    pof = mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma)
    toto = c(toto,pof)
  }
  E_A <- mean(toto)
  return(1 - alpha - E_A)
}

FeasCount <- function(g,n_g)
{
  G = matrix(g, ncol = n_g)
  Gfeas = G<=0.
  Allfeas = array(1,dim(G)[1])
  for(i in 1:(dim(G)[2]))
  {
    Allfeas = Allfeas * Gfeas[,i]
  }
  
  return (length(which(Allfeas == 1))/dim(G)[1])
}
###################################################################################### 
########################################### Probability of the probality of violation
###################################################################################### 
Probability_C <- function(x,modelConst,alea,beta,Nsim,d,m,n_g)
{

  input <- matrix(0,dim(alea)[1],d+m)
  input[,((d+1):(d+m))] <- alea
  for(j in 1:dim(alea)[1]) input[j,1:d] <- x 
  
  inp = data.frame()
  for (lv in levels(modelConst$X$output))
  {
    inp = rbind(inp,data.frame(input,rep(lv,dim(alea)[1])))
  }
  inp = setNames(inp, modelConst$inputNames)
  
  # sim = kergp::simulate(modelConst,nsim=Nsim,newdata=inp,cond=TRUE)$sim 
  # sim = simulate_fix(modelConst,nsim=Nsim,newdata=inp,cond=TRUE)$sim
  sim = simulate_new(modelConst,nsim=Nsim,newdata=inp,cond=TRUE, nugget.sim = 1e-10)
  
  sim = t(sim)
  # test <- apply(sim,1,function(vect) length(which(vect<=0.))/dim(inp)[1])
  test <- apply(sim,1,FeasCount,n_g)
  
  c <- length(which(test >= beta))/Nsim
  return(c=c)
  
}


###################################################################################### 
########################################### current feasible minimum
###################################################################################### 
minimin_g <- function(model,modelConstraint,d,m,alea,beta,n_g)
{
  X = model@X[,1:d]
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  res_m = future_apply(X,1,IntegratedProcess,model,alea,d,m,select="TRUE")
  
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  res_P = future_apply(X,1,Expectation_C,modelConstraint,alea,d,m,n_g)
  
  index <- which(res_P <= 0 )
  if(length(index) > 0 ) {index2 <- which.min(res_m[index]);   return(list(min=res_m[index[index2]],opt=index[index2]))}
  if(length(index) == 0 ) {t <- which.min(res_P); return(list(min=res_m[t],opt=t))}
}


###################################################################################### 
########################################### Feasible Expected Improvement 
###################################################################################### 
Feas_Expected_Improvement <- function(x,minimum_g,model,alea,d,m,Nsim,modelConst,beta,n_g)
{
  sdZ <- abs(IntegratedProcess(x,model,alea,d,m,select="FALSE")) #sqrt
  espZ <- (IntegratedProcess(x,model,alea,d,m,select="TRUE"))
  v <- (minimum_g-espZ)/sdZ
  phi <- dnorm(v, mean = 0, sd = 1, log = FALSE)
  PHI <- pnorm(v, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)  
  EIZ <- sdZ*(v*PHI+phi)
  proba <- Probability_C(x,modelConst,alea,beta,Nsim,d,m, n_g)
  return(c(EIZ)*proba ) # beware give negative values: -EIZ
}


###################################################################################### 
########################################### Variance of the integrated process at
########################################### "xstar" when we add "xp1" (step t+1)
###################################################################################### 
PredVar_MC_Z <- function(xstar,xp1,model,alea,d,m)
{
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
  X1 <- cbind(t(matrix(rep(xstar,nrow(alea)),length(xstar),nrow(alea))),alea)
  X2 <- (as.matrix( rbind(model@X,as.vector(xp1))))
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  dim <- ncol(X1)
  object <- objectinit
  object@range.val <- covparam2vect(objectinit)
  M <- covMat1Mat2(model@covariance,X1,X2,nugget.flag = TRUE)/model@covariance@sd2
  X1 <- (as.matrix( rbind(model@X,as.vector(xp1))))
  X2 <- cbind(t(matrix(rep(xstar,nrow(alea)),length(xstar),nrow(alea))),alea) 
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

###################################################################################### 
########################################### variance at "t+1" of the process in (XxU) space
########################################### at "newdata" when we add the point "xp1"
###################################################################################### 
PredVar_var_F_new <- function(newdata,xp1,model,d,m)
{
  n = dim(newdata)[1]+1
  utile <- predict(model,rbind(xp1,newdata),checkNames = FALSE,type="SK", covCompute = TRUE)
  return(utile$sd[2:(n)]^2 - utile$cov[2:(n),1]^2/utile$sd[1]^2)
}



###################################################################################### 
########################################### Mean at "t+1" of the integrated process in X-space
########################################### at "x" when we add the point "xp1" and "observation"
###################################################################################### 
Pred_mean_Z <- function(x,xp1,observation,model,alea,d,m)
{
  m1 <- matrix(x, nrow = 1, ncol = d)
  m2 = t(matrix(rep(x,dim(alea)[1]),d,dim(alea)[1]))
  dat <- data.frame(cbind(m2,alea))
  colnames(dat) <- NULL
  xp1 <- as.vector(xp1)
  data <- data.frame(rbind(xp1,dat))
  utile <- DiceKriging::predict.km(model,data,checkNames = FALSE,type="SK",cov.compute=TRUE)
  term1 <- IntegratedProcess(x,model,alea,d,m,"TRUE")
  term2 <- (observation - utile$mean[1])/utile$sd[1]^2
  term3 <- mean(utile$cov[1,])  
  return(term1 + term2*term3)
}  

###################################################################################### 
########################################### Variance of the Improvement at step "t+1"
########################################### when we add the point "xp1" and "observation"
###################################################################################### 
TotalVar <- function(x,xp1,model,alea,d,m,observation,ming,variZplus1)
{
  meanZplus1 <- Pred_mean_Z(x=x,xp1=xp1,observation,model=model,alea=alea,d=d,m=m)
  term <- (ming - meanZplus1)/sqrt(variZplus1)
  a <- pnorm(term, mean = 0, sd = 1)
  b <- dnorm(term, mean = 0, sd = 1)
  
  term1 <- (ming - meanZplus1)^2 + sqrt(variZplus1)^2
  term2 <- sqrt(variZplus1)*(ming - meanZplus1)
  
  return(c((ming - meanZplus1)*a + sqrt(variZplus1)*b,
           term1*a + term2*b - ((ming - meanZplus1)*a + sqrt(variZplus1)*b)^2 ))
}

###################################################################################### 
########################################### Total Variance of the improvement
########################################### at "xnew" when we add the point "xp1" and "observation"
###################################################################################### 
unext <- function(xp1,xnew,rep,modelObjective,alea,ming,m,d)
{
  oldlaw <- DiceKriging::predict.km(modelObjective,rbind(xp1),checkNames = FALSE,type="SK",cov.compute=TRUE)
  oldlaw <- c(oldlaw$mean,oldlaw$sd)
  aleaa <- rnorm(5000,oldlaw[1],max(oldlaw[2],1e-6))
  aleaa <- choice.grid(aleaa,rep,ng=1)
  aleaa <- as.vector(aleaa$opti_grid[,1])
  # aleaa <- matrix(randtoolbox::sobol(n=rep, dim = 1, init = TRUE, scrambling = 0, seed = 4711, normal = TRUE),ncol=1)
  # aleaa <- aleaa*oldlaw[2] + oldlaw[1]
  variZplus1 <- PredVar_MC_Z(xnew,xp1,modelObjective,alea,d,m)
  stopping <- 0
  term1 <- term2 <- NULL
  repeat{
    stopping = stopping+1
    muet <- TotalVar(xnew,xp1,modelObjective,alea,d,m,aleaa[stopping],ming,variZplus1)
    term1 <- c(term1,muet[1]) ;   term2 <- c(term2,muet[2])
    if (stopping == rep) break
  }
  return(var(term1) + mean(term2))
}

###################################################################################### 
########################################### Sampling
###################################################################################### 

Sampling_new_singleu <- function(xnew,rep,modelObjective,alea,ming,inputsss,modelConstraint,m,d,d_c)
{
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  
  res <- future_apply(inputsss,1,unext,xnew=xnew,rep=rep,modelObjective = modelObjective,alea = alea,ming = ming,m = m,d = d)
  
  n = dim(inputsss)[1]
  
  predvec_base = data.frame()
  for (cat in 1:d_c)
  {
    predvec_base = rbind(predvec_base, data.frame(inputsss,factor(array(levels(modelConstraint$X$output)[cat],dim(inputsss)[1]))))
  }
  predvec_base = setNames(predvec_base,modelConstraint$inputNames)
  
  pp = c()
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  pp = future_apply(inputsss,1,AuxFunSingleu,predvec_base,inputsss,modelConstraint,m,d,d_c)
  
  
  newxu = list()
  for(const in 1:d_c) newxu[[const]] = predvec_base[(n*(const-1)+which.min(res*pp)),]
  return(newxu)
}

Sampling_new_singleg <- function(xnew,rep,modelObjective,alea,ming,inputsss,modelConstraint,m,d,d_c)
{
  registerDoFuture()
  plan(multisession, workers = max(detectCores()-1, 1))
  
  res <- future_apply(inputsss,1,unext,xnew=xnew,rep=rep,modelObjective = modelObjective,alea = alea,ming = ming,m = m,d = d)
  
  n = dim(inputsss)[1]
  
  predvec_base = data.frame()
  for (cat in 1:d_c)
  {
    predvec_base = rbind(predvec_base, data.frame(inputsss,factor(array(levels(modelConstraint$X$output)[cat],dim(inputsss)[1]))))
  }
  predvec_base = setNames(predvec_base,modelConstraint$inputNames)
  
  pp = c()
  for(const in 1:d_c)
  {
    registerDoFuture()
    plan(multisession, workers = max(detectCores()-1, 1))
    ppl = future_apply(inputsss,1,AuxFunSingleg,predvec_base,inputsss,modelConstraint,m,d,d_c,const)
    pp = c(pp,ppl)
  }
  newxu = list()
  newxu[[1]] = inputsss[which.min(res),]
  newxu[[2]] = predvec_base[which.min(pp),]
  
  return(newxu)
}

AuxFunSingleg = function(input4s,predvec_base,inputsss,modelConstraint,m,d,d_c,const)
{
  
  input4smatrix = matrix(0,nrow = d_c, ncol = (d+m))
  for (cat in 1:d_c)
  {
    input4smatrix[cat,] = input4s
  }
  predvec_xp1 = data.frame(input4smatrix,levels(modelConstraint$X$output))
  predvec_xp1 = setNames(predvec_xp1,modelConstraint$inputNames)
  predvec = rbind(predvec_base,predvec_xp1 )
  
  
  utile <- predict(modelConstraint,predvec,checkNames = FALSE,type="SK",covCompute =TRUE)
  
  Knew = utile$cov[(d_c*n+const),(d_c*n+const)] 
  pq = c()
  for(i in 1:n)
  {
    pred_mean = c()
    for (cat in 1:d_c)
    {
      pred_mean = c(pred_mean,utile$mean[(n*(cat-1)+i)])
    }
    
    cov_l = matrix(0, nrow = d_c, ncol = d_c)
    for(j in 1:d_c)
    {
      for(k in 1:d_c)
      {
        
        kn = utile$cov[(i+(j-1)*n),(d_c*n+const)] 
        knp = utile$cov[(d_c*n+const),(i+(k-1)*n)] 
        
        cov_l[j,k] = utile$cov[(i+(j-1)*n),(i+(k-1)*n)]  - kn*knp/Knew
        cov_l = cov_l + diag(d_c)*1e-6
      }
    }
    
    lower <- rep(-Inf, d_c)
    upper <- rep(0, d_c)
    pof = mvtnorm::pmvnorm(lower = lower, upper = upper, mean = pred_mean, sigma = cov_l)
    
    pq = c(pq,pof)
  }
  
  pp <- mean(pq*(1-pq))
  return(pp)
}

AuxFunSingleu = function(input4s,predvec_base,inputsss,modelConstraint,m,d,d_c)
{
  
  input4smatrix = matrix(0,nrow = d_c, ncol = (d+m))
  for (cat in 1:d_c)
  {
    input4smatrix[cat,] = input4s
  }
  predvec_xp1 = data.frame(input4smatrix,levels(modelConstraint$X$output))
  predvec_xp1 = setNames(predvec_xp1,modelConstraint$inputNames)
  predvec = rbind(predvec_base,predvec_xp1 )
  
  
  utile <- predict(modelConstraint,predvec,checkNames = FALSE,type="SK",covCompute =TRUE)
  
  Knew = utile$cov[(d_c*n+1):(d_c*n+d_c),(d_c*n+1):(d_c*n+d_c)] 
  pq = c()
  for(i in 1:n)
  {
    pred_mean = c()
    for (cat in 1:d_c)
    {
      pred_mean = c(pred_mean,utile$mean[(n*(cat-1)+i)])
    }
    
    cov_l = matrix(0, nrow = d_c, ncol = d_c)
    for(j in 1:d_c)
    {
      for(k in 1:d_c)
      {
        kn = matrix(0, nrow = d_c, ncol = 1)
        knp = matrix(0, nrow = d_c, ncol = 1)
        
        for(ii in 1:d_c)
        {
          kn[ii] = utile$cov[(i+(j-1)*n),(d_c*n+ii)] 
          knp[ii] = utile$cov[(d_c*n+ii),(i+(k-1)*n)] 
        }
        cov_l[j,k] = utile$cov[(i+(j-1)*n),(i+(k-1)*n)]  - t(kn)%*%solve((Knew+diag(d_c)*1e-6),knp)
      }
    }
    
    lower <- rep(-Inf, d_c)
    upper <- rep(0, d_c)
    pof = mvtnorm::pmvnorm(lower = lower, upper = upper, mean = pred_mean, sigma = cov_l)
    
    pq = c(pq,pof)
  }
  
  pp <- mean(pq*(1-pq))
  
  return(pp)
}
