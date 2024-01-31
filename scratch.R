# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

simZ <- complement(my.data$Age, 0.8)
plot(my.data$Age, simZ)

cor(my.data$Age, simZ)

# Summary of results.===========================================================


lab.font.size <- 1.3
ax.font.size <- 1.5
title.font.size <- 1.8
# Box plot for difference in group sizes.
par(mfrow = c(1, 1))

pdf(file = "SpasmsCorr_DisImbalance.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(Balance~Method, data = all_metrics, 
        ylab = 'Absolute difference in group size',
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(Balance~Method, data = all_metrics, 
        ylab = 'Absolute difference in group size',
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
dev.off()

pdf(file = "SpasmsCorr_Energy.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(Energy~Method, data = all_metrics, 
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(Energy~Method, data = all_metrics, 
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
abline(h = original.result["Energy"], col = 'red',
       lty=2, lwd=3)
dev.off()

# Analysis of Age

pdf(file = "SpasmsCorr_Age_Means.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.Age~Method, data = all_metrics, 
        ylab = 'Difference Means Age', ylim = c(0, max(all_metrics$M.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(M.Age~Method, data = all_metrics, 
        ylab = 'Difference Means Age', ylim = c(0, max(all_metrics$M.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
abline(h = original.result["M.Age"], col = 'red',
       lty=2, lwd=3)
dev.off()

pdf(file = "SpasmsCorr_Age_SDs.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(SD.Age~Method, data = all_metrics, 
        ylab = 'Difference SD Age', ylim = c(0, max(all_metrics$SD.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(SD.Age~Method, data = all_metrics, 
        ylab = 'Difference SD Age', ylim = c(0, max(all_metrics$SD.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
abline(h = original.result["SD.Age"], col = 'red',
       lty=2, lwd=3)
dev.off()


# Analysis of ISFreq

pdf(file = "SpasmsCorr_ISFreq_Means.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.ISFreq~Method, data = all_metrics, 
        ylab = 'Difference Means FIS', ylim = c(0, max(all_metrics$M.ISFreq)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(M.ISFreq~Method, data = all_metrics, 
        ylab = 'Difference Means FIS', ylim = c(0, max(all_metrics$M.ISFreq)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
abline(h = original.result["M.ISFreq"], col = 'red',
       lty=2, lwd=3)
dev.off()

pdf(file = "SpasmsCorr_ISFreq_SDs.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(SD.ISFreq~Method, data = all_metrics, 
        ylab = 'Difference SD FIS', ylim = c(0, max(all_metrics$SD.ISFreq)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(SD.ISFreq~Method, data = all_metrics, 
        ylab = 'Difference SD FIS', ylim = c(0, max(all_metrics$SD.ISFreq)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
abline(h = original.result["SD.ISFreq"], col = 'red',
       lty=2, lwd=3)
dev.off()

# Analysis of correct guess probability

pdf(file = "SpasmsCorr_CG_Mean.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.CG~Method, data = all_metrics, 
        ylab = 'Mean Correct Guess Probability', ylim = c(0, max(all_metrics$M.CG)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(M.CG~Method, data = all_metrics, 
        ylab = 'Mean Correct Guess Probability', ylim = c(0, max(all_metrics$M.CG)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
abline(h = 0.5, col = 'red', lty=2, lwd=3)
dev.off()


# Function to simulate correlated variable.-------------------------------------

# The function called complement is taken from:
# https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

n.sam <- 100
energy.vec <- rep(NA, n.sam)
for (i in 1:n.sam){
  ind.A <- sample(1:n, n/2)
  ind.B <- !(1:n %in% ind.A)
  energy.vec[i] <- energy(Z[ind.A, ], Z[ind.B, ])
}
boxplot(energy.vec)

# ===================================================================
# Comparisons of three adaptive randomization methods for
# continuous covariates and two treatments.
#
# ===================================================================
library(e1071) # For "moments" function.

source("adaptive_randomization.R")
load("SimulatedDesignsCorrelated_Individual_V2.RData")


# Develop simulation protocol.
Z <- read.csv("data/Stand_Simulated_Corr.csv")
Z <- Z[,2:4]
p <- ncol(Z)
n <- nrow(Z)

# Set name for columns in the resulting compiled matrix.
my.col.names <- c("Method", "Balance", "M.Z1", "M.Z2", "M.Z3",
                  "SD.Z1", "SD.Z2", "SD.Z3", "Energy",
                  "M.CG", "SD.CG")

# COMPUTE PERFORMANCE METRICS.

# Ma and Hu (2013). Randomization algorithm based on Kernel Density.============
metrics <- lapply(results_kr, 
                  FUN = "performance.metrics.vec",
                  Z = Z)

metrics <- do.call(rbind, metrics)
metrics_kr <- data.frame("MH", metrics)
colnames(metrics_kr) <- my.col.names

# CA-RO Algorithm of Bertsimas et al. (2019).===================================
metrics <- lapply(results_caro, FUN = "performance.metrics.vec",
                  Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro <- data.frame("BKW", metrics)
colnames(metrics_caro) <- my.col.names

# Random allocation.============================================================
metrics <- lapply(results_random, FUN = "performance.metrics.vec",
                  Z = Z)
metrics <- do.call(rbind, metrics)
metrics_rand <- data.frame("RAND", metrics)
colnames(metrics_rand) <- my.col.names  

# Minimization algorithm for mean and variances of Nishi and Takaichi (2003).===
metrics <- lapply(results_mv, FUN = "performance.metrics.vec",
                  Z = Z)
metrics <- do.call(rbind, metrics)
metrics_mv <- data.frame("NT",metrics)
colnames(metrics_mv) <- my.col.names  

# Pocock and Simon method.======================================================
metrics <- lapply(results_ps, FUN = "performance.metrics.vec",
                  Z = Z)
metrics <- do.call(rbind, metrics)
metrics_ps <- data.frame("PS", metrics)
colnames(metrics_ps) <- my.col.names

# Save results.
all_metrics <- rbind(metrics_kr, metrics_caro, metrics_mv, metrics_ps)
#write.csv(all_metrics, file = "all_metrics_simulated_pembro.csv")

# Summary of results.===========================================================

lab.font.size <- 1.3
ax.font.size <- 1.5
title.font.size <- 1.8
# Box plot for difference in group sizes.
par(mfrow = c(1, 1))

#pdf(file = "SimCorr_DisImbalance.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(Balance~Method, data = all_metrics, 
        ylab = 'Absolute difference in group size',
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(Balance~Method, data = all_metrics, 
        ylab = 'Absolute difference in group size',
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
#dev.off()

#pdf(file = "SimCorr_Energy.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(Energy~Method, data = all_metrics, 
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(Energy~Method, data = all_metrics, 
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
abline(h = min(energy.dis.vec), col = 'red',lty=2, lwd=3)
#dev.off()


# Analysis of correct guess probability

#pdf(file = "SimCorr_CG_Mean.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.CG~Method, data = all_metrics, 
        ylab = 'Mean Correct Guess Probability', ylim = c(0, max(all_metrics$M.CG)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(M.CG~Method, data = all_metrics, 
        ylab = 'Mean Correct Guess Probability', ylim = c(0, max(all_metrics$M.CG)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
abline(h = 0.5, col = 'red', lty=2, lwd=3)
#dev.off()


# Compute optimal allocation for energy distance.===============================
# NOT TO RUN
library(snow)

n.comb.two <- combn(n, n/2)
n.choose.two <- ncol(n.comb.two)
my.fun <- function(x, Z, n){
  source("auxiliar_functions.R")
  W.one <- Z[x, ]
  W.two <- Z[!(1:n %in% x), ]
  energy(W.one, W.two)
}

cl <- makeSOCKcluster(c("localhost","localhost", "localhost", "localhost", "localhost", "localhost",
                        "localhost", "localhost"))
init.time <- proc.time()
energy.dis.vec <- parApply(cl, n.comb.two[,1:(n.choose.two/2)], MARGIN = 2, 
                           FUN = my.fun, Z = Z, n = n)
final.time <- proc.time() - init.time
stopCluster(cl)
save(energy.dis.vec, file = "EnergyDistance_Vectors_SimsCorrelated_Pembr.RData")
#min(energy.dis.vec)
#> final.time
#user  system elapsed 
#0.988   0.363 203.352 

# ===================================================================
# Comparisons of five adaptive randomization methods for three
# continuous covariates and two treatments.
# The simulations involve three dependent random variables based on the
# Pembro clinical trial.
#
# Construction of designs using different strategies and parameter settings.
# 
# ===================================================================

library(randomizeR)
library(e1071) # For "moments" function.

source("adaptive_randomization.R")

max.sim <- 1000

# Set the parameters for all algorithms.
n.zero <- 8 
bc.prob <- 0.8

# Set parameters for Bertsimas et al's approach.
rho.p <- 6
upper.gamma.bound <- 4

# Set parameters for Pocock and Simons.
n.categories <- 3

# Develop simulation protocol.
p <- 3 # Number of covariates.
my.data <- read.table("data/Pembrolizumab.txt", header = T)
my.data <- my.data[order(my.data$ID),]

n <- nrow(my.data) # Number of observations.

X <- matrix(NA, ncol = p, nrow = n)
# Simulate variable with high correlation.
set.seed(2523316)
X[,1] <- my.data[, "Age"]
X[,2] <- rnorm(n, mean = X[,1], sd = 5)
X[,3] <- rexp(n, 1/X[,2])

# Standardize the matrix of covariates.
Z <- scale(X)

# Save matrices for reproducibility.
write.csv(X, file = "data/Simulated_Corr.csv")
write.csv(Z, file = "data/Stand_Simulated_Corr.csv")


# Ma and Hu (2013). Randomization algorithm based on Kernel Density.============
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_kr <- apply(matrix(1:max.sim), MARGIN = 1, 
                    FUN = "kernel.randomization",
                    Z = Z, n.zero = n.zero, 
                    bc.prob = bc.prob)



# CA-RO Algorithm of Bertsimas et al. (2019).===================================
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_caro <- apply(matrix(1:max.sim), MARGIN = 1, 
                      FUN = "ca.ro",
                      Z = Z, rho.p = rho.p,
                      upper.gamma.bound = upper.gamma.bound,
                      n.zero = n.zero, 
                      bc.prob = 1) 


# Random allocation.============================================================
set.seed(1389874)
results_random <- apply(matrix(1:max.sim), MARGIN = 1, 
                        FUN = "random", Z = Z, n.zero = n.zero) 


# Minimization algorithm for mean and variances of Nishi and Takaichi (2003).===
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_mv <- apply(matrix(1:max.sim), MARGIN = 1, 
                    FUN = "minimization.mv",
                    Z = Z, n.zero = n.zero,
                    bc.prob = bc.prob) 



# Pocock and Simon method.======================================================
set.seed(1389874)
# Apply covariate adaptive allocation method under different setting. 
results_ps<- apply(matrix(1:max.sim), MARGIN = 1, 
                   FUN = "PS.minimization",
                   Z = Z, n.categories = n.categories,
                   n.zero = n.zero, 
                   bc.prob = bc.prob) 


# Save all designs.=============================================================
save(results_kr, results_caro, results_mv, results_ps, results_random,
     file = "SimulatedDesignsCorrelated_Individual.RData")


## REMOVE LATER

# Compute optimal allocation for energy distance.===============================
# NOT TO RUN
library(snow)

n.comb.two <- combn(n, n/2)
n.choose.two <- ncol(n.comb.two)
my.fun <- function(x, Z, n){
  source("auxiliar_functions.R")
  W.one <- Z[x, ]
  W.two <- Z[!(1:n %in% x), ]
  energy(W.one, W.two)
}

cl <- makeSOCKcluster(c("localhost","localhost", "localhost", "localhost", "localhost", "localhost",
                        "localhost", "localhost"))
init.time <- proc.time()
energy.dis.vec <- parApply(cl, n.comb.two[,1:(n.choose.two/2)], MARGIN = 2, 
                           FUN = my.fun, Z = Z, n = n)
final.time <- proc.time() - init.time
stopCluster(cl)
save(energy.dis.vec, file = "EnergyDistance_Vectors_Pembr.RData")
min(energy.dis.vec)
