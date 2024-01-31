# ===================================================================
# Adaptive randomization methods for continuous covariates
# and two treatment groups.
# 
# Alan R. Vazquez
# University of California, Los Angeles
# ===================================================================

library(randomizeR) # For "pbrPar" permuted randomized block.
library(Minirand) # For Pocock and Simon minimization method.
source("auxiliar_functions.R")

# Ma and Hu (2013). Randomization algorithm based on Kernel Density.

kernel.randomization <- function(x, Z, n.zero = 8, bc.prob = 0.8,
                                w = NULL){
  #library(randomizeR) # For "pbrPar" permuted randomized block.
  #source("auxiliar_functions.R")
  
  # Set parameters.
  #n.zero <- 8 # Initial number of patients (must be divisible by 2)
  p <- ncol(Z) # Number of covariates.
  if ( is.null(w) ){w <- rep(1/p, times = p)} # Same weights for the covariates.
  #bc.prob <- 0.8 # Biased coin probability.
  n <- nrow(Z)
  
  # Set objects to record correct guess probability.
  CGi <- rep(NA, times = n - n.zero)
  
  # Step 1. Allocate the first n.zero patients to two treatment groups.===========
  
  # 1.a. Select first subjects and define the initial group sizes.----------------
  ind.pat.in.trial <- 1:n.zero
  cov.pat.in.trial <- Z[ind.pat.in.trial,]
  n.one <- n.zero/2
  n.two <- n.zero/2
  
  # 1.b. Assign first set of subjects using permuted block randomization.---------
  params <- pbrPar(c(n.one,n.two), groups = c("1","2"))
  rs <- genSeq(params)
  BRD <- getRandList(rs)
  
  group.one <- which(BRD == "1")
  group.two <- which(BRD == "2")
  
  # Step 2. Assign the next patient.==============================================
  for (i in (n.zero+1):n){
    
    # 2.a. Standardize the covariates of the patients already in the trial.-------
    # Note: We do not need to normalize because the function assumes that all 
    # covariates are on the same scale. This is because our simulations 
    # assume that the whole covariate matrix has been obtained.
    
    # 2.b. Calculate imbalance measure.---------------------------------------------
    # Calculate imbalance measure.
    delta.d <- rep(NA, p)
    for (j in 1:p){
      f.one.z <- kernel.est(Z[i,j], cov.pat.in.trial[group.one,j])
      f.two.z <- kernel.est(Z[i,j], cov.pat.in.trial[group.two,j])
      delta.d[j] <- (n.one/n)*f.one.z - (n.two/n)*f.two.z
    }
    IMB <- sum(w*delta.d, na.rm = T) # New: Remove missing observations.
    
    # Step 3a. Determine assignment using a biased coin strategy.==================
    treat.assigned <- determine.assignment(IMB, bc.prob)
    
    # Step 3b. Calculate the correct guess probability.===========================
    CGi[i - n.zero] <- correct_guess_prob(n.one, n.two, treat.assigned)
    
    # Step 4. Update the groups and subjects in trial.============================
    if (treat.assigned == 1){
      # Augment group 1.
      n.one <- n.one + 1
      group.one <- c(group.one, i)
    } else { 
      # Augment group 2.
      n.two <- n.two + 1
      group.two <- c(group.two, i)
    }
    cov.pat.in.trial <- Z[1:i,] 
  }
  
  group.assignment <- list()
  group.assignment$G.one <- group.one
  group.assignment$G.two <- group.two
  group.assignment$CGi <- CGi # Store correct guess probability.
  
  return(group.assignment)
}

# CA-RO Algorithm of Bertsimas et al. (2019)
ca.ro <- function(x, Z, rho.p = 6, upper.gamma.bound = 4, n.zero = 8, bc.prob = 0.8){
  
  #library(randomizeR) # For "pbrPar" permuted randomized block.
  #source("auxiliar_functions.R")
  
  # Set parameters.
  n <- nrow(Z) # Number of subjects.
  p <- ncol(Z) # Number of covariates.
  gamma.bounds <- c(0.5, upper.gamma.bound)
  
  k <- n/2
  
  # Set objects to record correct guess probability.
  CGi <- rep(NA, times = n - n.zero)
  
  # Step 1. Allocate the first n.zero patients to two treatment groups.===========
  
  # 1.a. Select first subjects and define the initial group sizes.----------------
  n.one <- n.zero/2
  n.two <- n.zero/2
  ref.vector <- c(1, 2)
  
  # 1.b. Assign first set of subjects using permuted block randomization.---------
  params <- pbrPar(c(n.one,n.two), groups = c("1","2"))
  rs <- genSeq(params)
  BRD <- getRandList(rs)
  group.one <- which(BRD == "1")
  group.two <- which(BRD == "2")
  
  # Step 2. Calculate the discrepancy measure=====================================
  for (i in (n.zero+1):n){
    
    # 2.a. Check if it is feasible to assign treatments.
    bool.feasible.assignments <- c(n.one + 1 <= k, n.two + 1 <= k)
    
    # If it is feasible to add a subject to either group.
    if (all(bool.feasible.assignments)){ 
      
      # 2.a. Compute the absolute difference in means and variances of the two groups. 
      # Define uncertainty parameter at random.
      gamma.p <- runif(1, min = gamma.bounds[1], max = gamma.bounds[2])
      gamma.tilde <- (gamma.p^2)*(n - i)*p
      
      # 2.b. Compute the discrepancy measure for two groups.----------------------
      disc.groups <- bertsimas.discrepancy(i, Z, p, k, group.one, group.two, n,
                                           n.one, n.two, rho.p, gamma.tilde)
      
      # Step 2.c. Determine assignment using a biased coin strategy.==================
      treat.assigned <- determine.assignment(disc.groups[1] - disc.groups[2], bc.prob)
      
      
    } else { # If one group is full, then assign to the feasible group.
      treat.assigned <- ref.vector[bool.feasible.assignments]
    }
    
    # Step 2.d. Calculate the correct guess probability.===========================
    CGi[i  - n.zero] <- correct_guess_prob(n.one, n.two, treat.assigned)
    
    # Step 3. Update the groups and subjects in trial.============================
    if (treat.assigned == 1){
      # Augment group 1.
      n.one <- n.one + 1
      group.one <- c(group.one, i)
    } else { 
      # Augment group 2.
      n.two <- n.two + 1
      group.two <- c(group.two, i)
    }
    
  } # End for
  
  group.assignment <- list()
  group.assignment$G.one <- group.one
  group.assignment$G.two <- group.two
  group.assignment$CGi <- CGi # Store correct guess probability.
  
  return(group.assignment)
}

# Random allocation of patients.

random <- function(x, Z, n.zero = 8){
  
  #library(randomizeR) # For "pbrPar" permuted randomized block.
  #source("auxiliar_functions.R")
  n <- nrow(Z)
  
  # Set objects to record correct guess probability.
  CGi <- rep(NA, times = n - n.zero)
  
  # 1. Assign first set of subjects using permuted block randomization.---------
  n.one <- n.zero/2
  n.two <- n.zero/2
  
  params <- pbrPar(c(n.one,n.two), groups = c("1","2"))
  rs <- genSeq(params)
  BRD <- getRandList(rs)
  group.one <- which(BRD == "1")
  group.two <- which(BRD == "2")
  
  for (i in (n.zero+1):n){
    # Assign patient at random.
    treat.assigned <- sample(c(1,2), 1, replace = FALSE)
    
    # Step 3. Calculate the correct guess probability.===========================
    CGi[i  - n.zero] <- correct_guess_prob(n.one, n.two, treat.assigned)
    
    # Step 4. Update the groups and subjects in trial.============================
    if (treat.assigned == 1){
      # Augment group 1.
      n.one <- n.one + 1
      group.one <- c(group.one, i)
    } else { 
      # Augment group 2.
      n.two <- n.two + 1
      group.two <- c(group.two, i)
    }
  }
  
  group.assignment <- list()
  group.assignment$G.one <- group.one
  group.assignment$G.two <- group.two
  group.assignment$CGi <- CGi # Store correct guess probability.
  
  return(group.assignment)
}

# Minimization algorithm for mean and variances of Nishi and Takaichi (2003)

minimization.mv <- function(x, Z, probs = c(0.5, 0.5), n.zero = 8, bc.prob = 0.8){
  
  #library(randomizeR) # For "pbrPar" permuted randomized block.
  #source("auxiliar_functions.R")
  
  # Set parameters.
  p <- ncol(Z) # Number of covariates.
  n <- nrow(Z) # Number of subjects.
  
  # Set objects to record correct guess probability.
  CGi <- rep(NA, times = n - n.zero)
  
  # Step 1. Allocate the first n.zero patients to two treatment groups.===========
  
  # 1.a. Define the initial group sizes.----------------
  n.one <- n.zero/2
  n.two <- n.zero/2
  
  # 1.b. Assign first set of subjects using permuted block randomization.---------
  params <- pbrPar(c(n.one,n.two), groups = c("1","2"))
  rs <- genSeq(params)
  BRD <- getRandList(rs)
  group.one <- which(BRD == "1")
  group.two <- which(BRD == "2")
  
  # Step 2. Calculate the discrepancy measure=====================================
  for (i in (n.zero+1):n){
    
    # 2.b. Compute the discrepancy measure for two groups.----------------------
    disc.groups <- nishi.takaichi.discrepancy(i, Z, group.one, group.two, n.one, 
                                              n.two, probs)
    
    # Step 3a. Determine assignment using a biased coin strategy.==================
    treat.assigned <- determine.assignment(disc.groups[1] - disc.groups[2], bc.prob)
    
    # Step 3b. Calculate the correct guess probability.===========================
    CGi[i  - n.zero] <- correct_guess_prob(n.one, n.two, treat.assigned)
    
    # Step 4. Update the groups and subjects in trial.============================
    if (treat.assigned == 1){
      # Augment group 1.
      n.one <- n.one + 1
      group.one <- c(group.one, i)
    } else { 
      # Augment group 2.
      n.two <- n.two + 1
      group.two <- c(group.two, i)
    }
    
  } # End for
  
  group.assignment <- list()
  group.assignment$G.one <- group.one
  group.assignment$G.two <- group.two
  group.assignment$CGi <- CGi # Store correct guess probability.
  
  return(group.assignment)
}

# Pocock and Simon method.

PS.minimization <- function(x, Z, n.zero = 8, bc.prob = 0.8, n.categories = 3, w = NULL){
  
  #library(randomizeR) # For "pbrPar" permuted randomized block.
  #library(Minirand) # For Pocock and Simon minimization method.
  #source("auxiliar_functions.R")
  
  # Set parameters.
  p <- ncol(Z) # Number of covariates.
  n <- nrow(Z) # Number of subjects
  if ( is.null(w) ){w <- rep(1/p, times = p)} # Same weights for the covariates.
  
  # Set objects to record correct guess probability.
  CGi <- rep(NA, times = n - n.zero)
  
  # Step 0. Transform continuous covariates to discrete.
  Zd <- discrete.covariates(Z, n.categories)
  
  # Step 1. Allocate the first n.zero patients to two treatment groups.===========
  
  # 1.a. Select first subjects and define the initial group sizes.----------------
  n.one <- n.zero/2
  n.two <- n.zero/2
  
  # 1.b. Assign first set of subjects using permuted block randomization.---------
  params <- pbrPar(c(n.one,n.two), groups = c("1","2"))
  rs <- genSeq(params)
  BRD <- getRandList(rs)
  
  res <- rep(NA, times = n)
  res[1:n.zero] <- as.numeric(BRD)
  
  # Step 2. Assign the next patient.==============================================
  for (i in (n.zero+1):n){
    
    # Compute the current sizes of the groups.
    n.one <- sum(res == 1, na.rm = T)
    n.two <- sum(res == 2, na.rm = T)
    
    res[i] <- Minirand(covmat = Zd, i, covwt = w, ratio = c(1,1),
                       ntrt = 2, trtseq = c(1,2), method = 'Range',
                       result = res, p = bc.prob)
    
    # Save the correct guess probability.
    CGi[i  - n.zero] <- correct_guess_prob(n.one, n.two, res[i])
  }
  
  group.assignment <- list()
  group.assignment$G.one <- which(res == 1)
  group.assignment$G.two <- which(res == 2)
  group.assignment$CGi <- CGi # Store correct guess probability.
  
  return(group.assignment)
}

