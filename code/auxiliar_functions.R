# ===================================================================
# Auxiliary functions for adaptive randomization methods for 
# continuous covariates and two treatment groups.
# 
# Alan R. Vazquez
# University of California, Los Angeles
# ===================================================================
library(e1071) # For "moments" function.
library(expm) # For "sqrtm" function
library(gtools) # For "quantcut" function. 
library(energy) # For "eqdist.etest" function.

# Scaling functions.------------------------------------------------------------
scale.bert <- function(X){
  # Input: X an n x p matrix of p covariates on n subjects.
  # Output: std.X an n x p standardized matrix of p covariates on n subjects.
  #         The standarization is using Bertsimas et al.
  
  est.mean <- colMeans(X)
  est.var <- cov(X)
  trans.cov <- solve(est.var)
  trans.cov.sqrt <- sqrtm(trans.cov)
  std.X <- matrix(NA, ncol = ncol(X), nrow = nrow(X))
  for (i in 1:nrow(X)){
    std.X[i,] <- trans.cov.sqrt%*%(X[i,] - est.mean)  
  }
  result <- list(std.X, est.mean, trans.cov.sqrt)
  return(result)
}

# Auxiliary functions for kernel.randomization.----------------------------------
gaussian.kernel <- function(x){( (2*pi)^{-1/2} )*exp( -1*(x^2)/2 )}
h.bandwith <- function(x, s){s*(x^{-0.2})}
kernel.est <- function(z, Z){
  # Inputs
  #   z: covariate value of new subject.
  #   Z: covariate vector of subjects in trial.
  #   
  n.k <- length(Z)
  h.val <- h.bandwith(n.k, sd(Z))
  kernel.sum <- sum( sapply((z - Z)/h.val, FUN = "gaussian.kernel") )
  return(kernel.sum/(n.k*h.val))
}
determine.assignment <- function(IMB, bc.prob = 0.8){
  if (IMB < 0){
    treat.assigned <- sample(c(1,2), 1, prob = c(bc.prob, 1 - bc.prob), replace = F)
  } else if (IMB > 0){
    treat.assigned <- sample(c(1,2), 1, prob = c(1 - bc.prob, bc.prob), replace = F)
  } else {
    treat.assigned <- sample(c(1,2), 1, prob = c(0.5, 0.5), replace = F)
  }
  return(treat.assigned)
}

# Auxiliary functions for PS.minimization.--------------------------------------

discrete.covariates <- function(Z, n.categories = 3){
  p <- ncol(Z)
  n <- nrow(Z)
  
  Zd <- matrix(NA, ncol = p, nrow = n)
  for (j in 1:p){
    Zd[,j] <- quantcut(Z[,j], q = n.categories)
  }
  
  return(Zd)
}

Pocock.Simon.discrepancy <- function(i, n.one){
  D.t <- 2*n.one - (i - 1)
  return(D.t)
}


# Auxiliary functions for ca.ro algorithm.--------------------------------------
Phi <- function(k, l, x.t){
  if (k - l - x.t >= 1){return(1)} else {return(0)}
}

bertsimas.discrepancy <- function(i, Z, p, k, group.one, group.two, n, n.one, n.two, 
                                  rho.p, gamma.tilde){
  w.t <- Z[i,] # Select covariate of the subject to be assigned.
  
  # Compute 
  W.bar <- colMeans(Z[1:i,])
  Sigma.cov <- cov(Z[1:i,])
  V <- sqrtm(Sigma.cov)
  
  # Difference in means.
  abs.dif.cov.means <- matrix(NA, nrow = p, ncol = 2)
  # Difference in variances.
  abs.dif.cov.varns <- matrix(NA, nrow = p, ncol = 2)
  
  for (j in 1:p){
    
    # Difference in means.
    first.term.m <- sum(Z[group.one,j] - W.bar[j]) -  sum(Z[group.two,j] - W.bar[j])
    third.term.m <- sqrt(gamma.tilde*(n-i)*sum(V[j,]^2))
    
    #abs.dif.cov.means[j,1] <- max(first.term.m + (w.t[j] - W.bar[j]) + third.term.m,
    #                              -1*first.term.m - (w.t[j] - W.bar[j]) + third.term.m)
    abs.dif.cov.means[j,1] <- abs(first.term.m + (w.t[j] - W.bar[j])) + third.term.m
    
    #abs.dif.cov.means[j,2] <- max(first.term.m - (w.t[j] - W.bar[j]) + third.term.m,
    #                              -1*first.term.m + (w.t[j] - W.bar[j]) + third.term.m)
    abs.dif.cov.means[j,2] <- abs(first.term.m - (w.t[j] - W.bar[j])) + third.term.m
    
    # Difference in variances.
    first.term.v <- sum( (Z[group.one,j] - W.bar[j])^2 ) - sum( (Z[group.two,j] - W.bar[j])^2 )
    third.term.v <- gamma.tilde*sum(V[j,]^2)
    
    abs.dif.cov.varns[j,1] <- max(first.term.v + (w.t[j] - W.bar[j])^2 + third.term.v*Phi(k, n.one, 1),
                                  -1*first.term.v - (w.t[j] - W.bar[j])^2 + third.term.v*Phi(k, n.two, 0))
    
    abs.dif.cov.varns[j,2] <- max(first.term.v - (w.t[j] - W.bar[j])^2 + third.term.v*Phi(k, n.one, 0),
                                  -1*first.term.v + (w.t[j] - W.bar[j])^2 + third.term.v*Phi(k, n.two, 1))
    
  }
  
  # Compute the discrepancy measure for two groups.
  # Create discrepancy matrix.
  # We use the sqrt of abs.dif.cov.varns so that the quantities are on the same
  # scale.
  disc.mat <- abs.dif.cov.means/k + rho.p*sqrt(abs.dif.cov.varns/k)
  disc.groups <- colSums(disc.mat)
  
  return(disc.groups)
}

# Auxiliary functions for minimization.mv.--------------------------------------

nishi.takaichi.discrepancy <- function(i, Z, group.one, group.two, n.one, n.two,
                                       probs){
  
  # Calculate group averages and standard deviations of covariates for each group.
  M.g.one <- colMeans(Z[group.one,])
  M.g.two <- colMeans(Z[group.two,])
  SD.g.one <- apply(Z[group.one,], 2, sd)
  SD.g.two <- apply(Z[group.two,], 2, sd)
  
  # Calculate total average and standard deviations of covariates across groups.
  N.t <- n.one + n.two
  M.t <- (n.one*M.g.one + n.two*M.g.two)/N.t # Weighted average.
  Var.t <- ( (n.one - 1)*(SD.g.one^2) + (n.two - 1)*(SD.g.two^2))/(N.t -2)
  SD.t <- sqrt(Var.t) # Weighted standard deviation.
  
  # Calculate revised group group averages and standard deviations of covariates
  # for each group. That is, the quantities that result from allocating the new
  # subject to a group.
  r.M.g.one <- colMeans(Z[ c(group.one, i),])
  r.M.g.two <- colMeans(Z[ c(group.two, i),])
  r.SD.g.one <- apply(Z[ c(group.one, i),], 2, sd)
  r.SD.g.two <- apply(Z[ c(group.two, i),], 2, sd)
  
  # Calculate revised total average and standard deviations of covariates across 
  # groups. That is, the quantities that result from including the new subject
  # in the trial.
  r.N.t <- n.one + n.two + 1
  
  # If we add subject to group 1.
  r.M.t.g.one <- ( (n.one+1)*r.M.g.one + n.two*M.g.two)/r.N.t # Weighted average.
  r.Var.t.g.one <- ( n.one*(r.SD.g.one^2) + (n.two - 1)*(SD.g.two^2))/(r.N.t - 2)
  r.SD.t.g.one <- sqrt(r.Var.t.g.one) # Weighted standard deviation.
  
  # If we add subject to group 2.
  r.M.t.g.two <- ( n.one*M.g.one + (n.two+1)*r.M.g.two)/r.N.t # Weighted average.
  r.Var.t.g.two <- ( (n.one - 1)*(SD.g.one^2) + n.two*(r.SD.g.two^2))/(r.N.t - 2)
  r.SD.t.g.two <- sqrt(r.Var.t.g.two) # Weighted standard deviation.
  
  # Calculate differences between M(g) and M.t and between r.M(g) and r.M.t
  disc.one.M <- abs(r.M.g.one - r.M.t.g.one) - abs(M.g.one - M.t)
  disc.two.M <- abs(r.M.g.two - r.M.t.g.two) - abs(M.g.two - M.t)
  
  # Calculate differences between SD(g) and SD.t and between r.SD(g) and r.SD.t
  disc.one.SD <- abs(r.SD.g.one - r.SD.t.g.one) - abs(SD.g.one - SD.t)
  disc.two.SD <- abs(r.SD.g.two - r.SD.t.g.two) - abs(SD.g.two - SD.t)
  
  # Calculate differences between proportion of patients.
  disc.one.N <- n.one/N.t - probs[1]
  disc.two.N <- n.two/N.t - probs[2]
  
  # Calculate overall discrepancy value for all covariates.
  # Warning: Equal weights to all measures.
  CS.one <- sum(disc.one.M + disc.one.SD) + disc.one.N
  CS.two <- sum(disc.two.M + disc.two.SD) + disc.two.N
  return(c(CS.one, CS.two))
}

# Performance metrics.----------------------------------------------------------

num.subject.diff <- function(W.one, W.two){
  n.one <- nrow(W.one)
  n.two <- nrow(W.two)
  return(abs(n.one - n.two))
}

means.dist <- function(W.one, W.two){
  m.dist <- norm(colMeans(W.one) - colMeans(W.two), type = "2")
  return(m.dist)
}

means.diff <- function(W.one, W.two){
  m.diff <- abs(colMeans(W.one) - colMeans(W.two))
  return(m.diff)
}

sd.dist <- function(W.one, W.two){
  s.dist <- norm(apply(W.one, 2, sd) - apply(W.two, 2, sd), type = "2")
  return(s.dist)
}

sd.diff <- function(W.one, W.two){
  s.diff <- abs(apply(W.one, 2, sd) - apply(W.two, 2, sd))
  return(s.diff)
}


energy <- function(W.one, W.two){
  
  # Apply twinning
  n.one <- nrow(W.one)
  n.two <- nrow(W.two)
  
  n.one.choose.two <- combn(1:n.one, m = 2)
  n.two.choose.two <- combn(1:n.two, m = 2)
  
  # Calculate ED.one
  ED.one <- matrix(NA, ncol = n.one, nrow = n.two)
  for (j in 1:n.two){
    ED.one[j,] <- sapply(1:n.one, FUN = function(x, U, v){
      norm(U[x,] - v, type = "2")}, U = W.one, v = W.two[j,])
  }
  
  # Calculate ED.two/2
  ED.two <- apply(n.one.choose.two, MARGIN = 2, 
                  FUN = function(x, U) {norm(U[x[1],] - U[x[2],], type = "2")}, 
                  U = W.one)
  
  # Calculate ED.three/2
  ED.three <- apply(n.two.choose.two, MARGIN = 2, 
                    FUN = function(x, U) {norm(U[x[1],] - U[x[2],], type = "2")}, 
                    U = W.two)
  
  e <- 2*sum(ED.one)/(n.one*n.two) - 2*sum(ED.two)/(n.one^2) - 2*sum(ED.three)/(n.two^2)  
  return(e)
}


performance.metrics <- function(X, Z){
  W.one <- Z[X$G.one,]
  W.two <- Z[X$G.two,]
  
  result <- list()
  # Joint distribution metrics
  result$n.subj.balance <- num.subject.diff(W.one, W.two)
  result$mean.difference <- means.diff(W.one, W.two)
  result$sd.difference <- sd.diff(W.one, W.two)
  result$energy <- energy(W.one, W.two)
  
  return(result)
}

correct_guess_prob <- function(n.one, n.two, treatment){
  if ( (n.one < n.two & treatment == 1) | (n.one > n.two & treatment == 2) ){CG <- 1}
  if (n.one == n.two){CG <- 0.5}
  if ( (n.one < n.two & treatment == 2) | (n.one > n.two & treatment == 1) ){CG <- 0}
  return(CG)
}

performance.metrics.vec <- function(X, Z){
  
  W.one <- Z[X$G.one,]
  W.two <- Z[X$G.two,]
  
  # Joint distribution metrics
  n.subj.balance <- num.subject.diff(W.one, W.two)
  mean.difference <- means.diff(W.one, W.two)
  sd.difference <- sd.diff(W.one, W.two)
  energy <- energy(W.one, W.two)
  mean_CGp <- mean(X$CGi, na.rm = TRUE)
  sd_CGp <- sd(X$CGi, na.rm = TRUE)
  
  result <- c(n.subj.balance, mean.difference, sd.difference, energy, mean_CGp, sd_CGp)
  return(result)
}

