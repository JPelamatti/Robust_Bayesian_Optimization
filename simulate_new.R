#Simulate trajectories of the Gaussian process

simulate_new <- function(object, nsim = 1, seed = NULL, newdata = NULL,
                        cond = TRUE, nugget.sim = 0, checkNames = TRUE, ...) {
  
  if (!is.numeric(nugget.sim)) stop("'nugget.sim' must be a number")
  if (nugget.sim<0) stop("nugget.sim (homogenous to a variance) must not be negative")
  if (!is.logical(cond)) stop("'cond' must be TRUE/FALSE")

  

    m <- dim(newdata)[1]

    if (!identical(colnames(newdata), colnames(object$X))) {
      colnames(newdata) <- colnames(object$X)
    }
    nms <- object$inputNames   
    New <- newdata[ , nms, drop = FALSE]
    tt <- delete.response(terms(object))
    mf <- model.frame(tt, data = data.frame(newdata))
    F.newdata <- model.matrix(tt, data = mf)
    XNew <- newdata[ , nms, drop = FALSE]
    Sigma <- covMat(object$covariance, X = XNew, Xnew = XNew,
                   compGrad = FALSE)
    T.newdata <- chol(Sigma + diag(nugget.sim, m, m))
  
  
    y.trend <- F.newdata %*% object$betaHat
  
    # simulations conditional to the observations
    #if (object@noise.flag) {
    #	stop("conditional simulations not available for heterogeneous observations")
    #} else {
    Sigma21 <- covMat(object$covariance, X = object$X, Xnew = XNew, compGrad = FALSE)          ## size n x m
    Tinv.Sigma21 <- backsolve(object$L, Sigma21, upper.tri = FALSE)     ## t(T22)^(-1) * Sigma21,  size  n x m
    y.trend.cond <- y.trend + t(Tinv.Sigma21) %*% object$eStar                 ## size m x 1
    
    Sigma11 <- Sigma	
    
    Sigma.cond <- Sigma11 - t(Tinv.Sigma21) %*% Tinv.Sigma21          ## size m x m
    T.cond <- chol(Sigma.cond + diag(nugget.sim, m, m))			
    white.noise <- matrix(rnorm(m*nsim), m, nsim)
    y.rand.cond <- t(T.cond) %*% white.noise
    y <- matrix(y.trend.cond, m, nsim) + y.rand.cond	

  
  return(y)
  
}