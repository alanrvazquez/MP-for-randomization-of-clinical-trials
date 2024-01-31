# ===================================================================
# Comparisons of three adaptive randomization methods for
# continuous covariates and two treatments.
#
# ===================================================================
library(e1071) # For "moments" function.

source("code/adaptive_randomization.R")
load("Study_Tuning_Parameters/PembroDesigns_CARO_Different_rho.RData")

# Develop simulation protocol.
p <- 3 # Number of covariates.
my.data <- read.table("data/Pembrolizumab.txt", header = T)
my.data <- my.data[order(my.data$ID),]

n <- nrow(my.data) # Number of observations.
# Standardize the matrix of covariates.
Z <- scale(data.matrix(my.data[,c("Age", "PD.L1.MPS", "Hemoglobin")]))

# Set name for columns in the resulting compiled matrix.
my.col.names <- c("Rho", "Balance", "M.Age", "M.PD.L1.MPS", "M.Hemoglobin",
                  "SD.Age", "SD.PD.L1.MPS", "SD.Hemoglobin", "Energy",
                  "M.CG", "SD.CG")

# COMPUTE PERFORMANCE METRICS.

# CA-RO Algorithm of Bertsimas et al. (2019).===================================
metrics <- lapply(results_caro_rhoOne, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_rhoOne <- data.frame("One", metrics)
colnames(metrics_caro_rhoOne) <- my.col.names

metrics <- lapply(results_caro_rhoThree, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_rhoThree <- data.frame("Three", metrics)
colnames(metrics_caro_rhoThree) <- my.col.names

metrics <- lapply(results_caro_rhoSix, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_rhoSix <- data.frame("Six", metrics)
colnames(metrics_caro_rhoSix) <- my.col.names

metrics <- lapply(results_caro_rhoTen, FUN = "performance.metrics.vec", Z = Z)
metrics <- do.call(rbind, metrics)
metrics_caro_rhoTen <- data.frame("Ten", metrics)
colnames(metrics_caro_rhoTen) <- my.col.names


# Save results.
all_metrics <- rbind(metrics_caro_rhoOne, metrics_caro_rhoThree,
                     metrics_caro_rhoSix, metrics_caro_rhoTen)
all_metrics$Rho <- factor(all_metrics$Rho, levels = c("One", "Three", "Six", "Ten"))
write.csv(all_metrics, file = "metrics_pembro_rho_CARO.csv")

all_metrics <- read.csv("metrics_pembro_rho_CARO.csv")
all_metrics$Rho <- factor(all_metrics$Rho, levels = c("One", "Three", "Six", "Ten"))

