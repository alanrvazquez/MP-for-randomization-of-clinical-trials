# ===================================================================
# Comparisons of three adaptive randomization methods for
# continuous covariates and two treatments.
#
# ===================================================================
library(e1071) # For "moments" function.

source("adaptive_randomization.R")
load("PembroDesigns_CARO_Different_rho.RData")

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


# Summary of results.===========================================================
lab.font.size <- 1.3
ax.font.size <- 1.5
title.font.size <- 1.8
# Box plot for difference in group sizes.

#pdf(file = "Pembro_Energy.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
par(mfrow = c(1, 1))
boxplot(Energy~Rho, data = all_metrics, 
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
#dev.off()

# Analysis of Age

#pdf(file = "Pembro_Age_Means.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
par(mfrow = c(1,2))
boxplot(M.Age~Rho, data = all_metrics, 
        ylab = 'Difference Means Age', ylim = c(0, max(all_metrics$M.Age, all_metrics$SD.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(SD.Age~Rho, data = all_metrics, 
        ylab = 'Difference SD Age', ylim = c(0, max(all_metrics$M.Age, all_metrics$SD.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
#dev.off()


# Analysis of PD.L1.MPS

#pdf(file = "Pembro_PD_Means.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
par(mfrow = c(1,2))
boxplot(M.PD.L1.MPS~Rho, data = all_metrics, 
        ylab = 'Difference Means PD-L1', ylim = c(0, max(all_metrics$M.PD.L1.MPS, all_metrics$SD.PD.L1.MPS)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)

#pdf(file = "Pembro_PD_SDs.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(SD.PD.L1.MPS~Rho, data = all_metrics, 
        ylab = 'Difference SD PD-L1', ylim = c(0, max(all_metrics$M.PD.L1.MPS, all_metrics$SD.PD.L1.MPS)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
#dev.off()

# Analysis of Hemoglobin

par(mfrow = c(1,2))
#pdf(file = "Pembro_Hemoglobin_Means.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.Hemoglobin~Rho, data = all_metrics, 
        ylab = 'Difference Means Hemoglobin', ylim = c(0, max(all_metrics$M.Hemoglobin, all_metrics$SD.Hemoglobin)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
#dev.off()

#pdf(file = "Pembro_Hemoglobin_SDs.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(SD.Hemoglobin~Rho, data = all_metrics, 
        ylab = 'Difference SD Hemoglobin', ylim = c(0, max(all_metrics$M.Hemoglobin, all_metrics$SD.Hemoglobin)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
#dev.off()

par(mfrow = c(1,2))
#pdf(file = "Pembro_CG_Mean.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.CG~Rho, data = all_metrics, 
        ylab = 'Mean Correct Guess Probability', ylim = c(0, max(all_metrics$M.CG, all_metrics$SD.CG)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
#dev.off()

#pdf(file = "Pembro_CG_SD.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(SD.CG~Rho, data = all_metrics, 
        ylab = 'Standard Deviation Correct Guess Probability', ylim = c(0, max(all_metrics$M.CG, all_metrics$SD.CG)),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
#dev.off()
