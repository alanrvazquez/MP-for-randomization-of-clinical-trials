# ===================================================================
# Comparisons of three adaptive randomization methods for
# continuous covariates and two treatments.
#
# ===================================================================
library(e1071) # For "moments" function.

source("code/adaptive_randomization.R")
load("Study_Tuning_Parameters/PembroDesigns_CARO_Different_Gamma.RData")

# Develop simulation protocol.
p <- 3 # Number of covariates.
my.data <- read.table("data/Pembrolizumab.txt", header = T)
my.data <- my.data[order(my.data$ID),]

n <- nrow(my.data) # Number of observations.
# Standardize the matrix of covariates.
Z <- scale(data.matrix(my.data[,c("Age", "PD.L1.MPS", "Hemoglobin")]))

# Set name for columns in the resulting compiled matrix.
my.col.names <- c("Gamma", "Balance", "M.Age", "M.PD.L1.MPS", "M.Hemoglobin",
                  "SD.Age", "SD.PD.L1.MPS", "SD.Hemoglobin", "Energy",
                  "M.CG", "SD.CG")

# COMPUTE PERFORMANCE METRICS.

# CA-RO Algorithm of Bertsimas et al. (2019).===================================
metrics <- lapply(results_caro_GOne, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_GOne <- data.frame("One", metrics)
colnames(metrics_caro_GOne) <- my.col.names

metrics <- lapply(results_caro_GTwo, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_GTwo <- data.frame("Two", metrics)
colnames(metrics_caro_GTwo) <- my.col.names

metrics <- lapply(results_caro_GFour, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_GFour <- data.frame("Four", metrics)
colnames(metrics_caro_GFour) <- my.col.names

metrics <- lapply(results_caro_GEight, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_GEight <- data.frame("Eight", metrics)
colnames(metrics_caro_GEight) <- my.col.names


# Save results.
all_metrics <- rbind(metrics_caro_GOne, metrics_caro_GTwo,
                     metrics_caro_GFour, metrics_caro_GEight)
all_metrics$Gamma <- factor(all_metrics$Gamma, levels = c("One", "Two", "Four", "Eight"))
write.csv(all_metrics, file = "Study_Tuning_Parameters/metrics_pembro_Gamma_CARO.csv")


###############################################
### SPASMS CLINICAL TRIAL
###############################################

source("code/adaptive_randomization.R")
load("Study_Tuning_Parameters/SpasmsDesigns_CARO_Different_Gamma.RData")

p <- 2 # Number of covariates.
my.data <- read.table("data/Infant_Spasms.txt", header = T)
my.data <- my.data[order(my.data$ID),]

n <- nrow(my.data) # Number of observations.
# Standardize the matrix of covariates.
Z <- scale(data.matrix(my.data[,c("Age", "ISFreq")]))

# Set name for columns in the resulting compiled matrix.
my.col.names <- c("Gamma", "Balance", "M.Age", "M.ISFreq", "SD.Age", "SD.ISFreq", 
                  "Energy", "M.CG", "SD.CG")

# COMPUTE PERFORMANCE METRICS.

# CA-RO Algorithm of Bertsimas et al. (2019).===================================
metrics <- lapply(results_caro_GOne, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_GOne <- data.frame("One", metrics)
colnames(metrics_caro_GOne) <- my.col.names

metrics <- lapply(results_caro_GTwo, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_GTwo <- data.frame("Two", metrics)
colnames(metrics_caro_GTwo) <- my.col.names

metrics <- lapply(results_caro_GFour, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_GFour <- data.frame("Four", metrics)
colnames(metrics_caro_GFour) <- my.col.names

metrics <- lapply(results_caro_GEight, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_GEight <- data.frame("Eight", metrics)
colnames(metrics_caro_GEight) <- my.col.names


# Save results.
all_metrics <- rbind(metrics_caro_GOne, metrics_caro_GTwo,
                     metrics_caro_GFour, metrics_caro_GEight)
all_metrics$Gamma <- factor(all_metrics$Gamma, levels = c("One", "Two", "Four", "Eight"))
write.csv(all_metrics, file = "Study_Tuning_Parameters/metrics_spasms_Gamma_CARO.csv")



