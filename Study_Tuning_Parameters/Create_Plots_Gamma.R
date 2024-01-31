source("code/adaptive_randomization.R")

all_metrics <- read.csv("Study_Tuning_Parameters/metrics_pembro_Gamma_CARO.csv")
all_metrics <- all_metrics[,-1]

all_metrics$Gamma <- factor(all_metrics$Gamma, levels = c("One", "Two", "Four", "Eight"))


# Summary of results.===========================================================
lab.font.size <- 1.3
ax.font.size <- 1.5
title.font.size <- 1.8
# Box plot for difference in group sizes.

pdf(file = "Study_Tuning_Parameters/Pembro_Gamma.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
par(mfrow = c(1, 2))
boxplot(Energy~Gamma, data = all_metrics, xlab = "Maximum Gamma Value",
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(Energy~Gamma, data = all_metrics, xlab = "Maximum Gamma Value",
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)

boxplot(M.CG~Gamma, data = all_metrics, xlab = "Maximum Gamma Value",
        ylab = 'Mean Correct Guess Probability', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(M.CG~Gamma, data = all_metrics, xlab = "Maximum Gamma Value",
        ylab = 'Mean Correct Guess Probability', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
dev.off()

## SPASMS CLINICAL TRIAL.=======================================================

all_metrics <- read.csv("Study_Tuning_Parameters/metrics_spasms_Gamma_CARO.csv")
all_metrics <- all_metrics[,-1]

all_metrics$Gamma <- factor(all_metrics$Gamma, levels = c("One", "Two", "Four", "Eight"))


# Summary of results.===========================================================
lab.font.size <- 1.3
ax.font.size <- 1.5
title.font.size <- 1.8
# Box plot for difference in group sizes.

pdf(file = "Study_Tuning_Parameters/Spasms_Gamma.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
par(mfrow = c(1, 2))
boxplot(Energy~Gamma, data = all_metrics, xlab = "Maximum Gamma Value",
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(Energy~Gamma, data = all_metrics, xlab = "Maximum Gamma Value",
        ylab = 'Energy', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)

boxplot(M.CG~Gamma, data = all_metrics, xlab = "Maximum Gamma Value",
        ylab = 'Mean Correct Guess Probability', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NULL, lty = 6)
boxplot(M.CG~Gamma, data = all_metrics, xlab = "Maximum Gamma Value",
        ylab = 'Mean Correct Guess Probability', ylim = c(0, 1),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE)
dev.off()

