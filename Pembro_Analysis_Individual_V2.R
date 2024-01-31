# ===================================================================
# Comparisons of three adaptive randomization methods for
# continuous covariates and two treatments.
#
# ===================================================================
library(e1071) # For "moments" function.
library(RColorBrewer) # For colorful plots.
library(fmsb) # For radar chart.

source("code/adaptive_randomization.R")
load("designs/PembroDesigns_Individual_v2.RData")


# Develop simulation protocol.
p <- 3 # Number of covariates.
my.data <- read.table("data/Pembrolizumab.txt", header = T)
my.data <- my.data[order(my.data$ID),]

n <- nrow(my.data) # Number of observations.
# Standardize the matrix of covariates.
Z <- scale(data.matrix(my.data[,c("Age", "PD.L1.MPS", "Hemoglobin")]))

# Set name for columns in the resulting compiled matrix.
my.col.names <- c("Method", "Balance", "M.Age", "M.PD.L1.MPS", "M.Hemoglobin",
                  "SD.Age", "SD.PD.L1.MPS", "SD.Hemoglobin", "Energy",
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
all_metrics <- rbind(metrics_kr, metrics_caro,  metrics_mv, metrics_ps)


# Evaluate the actual assignment.===============================================

original.allocation <- list()
original.allocation$G.one <- which(my.data$Arm == 'A')
original.allocation$G.two <- which(my.data$Arm == 'B')
original.allocation$CGi <- rep(0, n)
original.result <- performance.metrics.vec(original.allocation, Z)
names(original.result) <- my.col.names[-1]

# Summary of results.===========================================================

# Set font sizes.
lab.font.size <- 1.3
ax.font.size <- 1.5
title.font.size <- 1.8

# Define color palette.
coul <- brewer.pal(5, "Set1")
colors_border <- coul
colors_in <- alpha(coul,0.3)
ref_line_col <- "black" # Color of reference line.

# Box plot for difference in group sizes.
par(mfrow = c(1, 1))

pdf(file = "figures/Pembro_DisImbalance.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(Balance~Method, data = all_metrics, 
        ylab = 'Absolute difference in group size', col = colors_in, 
        medcol = colors_border, outcol = colors_in,
        cex.lab = lab.font.size, cex.axis = ax.font.size, notch = TRUE)
grid(NULL,NULL, lty = 6)
boxplot(Balance~Method, data = all_metrics, 
        ylab = 'Absolute difference in group size', col = colors_in, 
        medcol = colors_border, outcol = colors_in,
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, notch = TRUE)
dev.off()

pdf(file = "figures/Pembro_Energy.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(Energy~Method, data = all_metrics, 
        ylab = 'Energy', ylim = c(0, 1), col = colors_in, 
        medcol = colors_border, outcol = colors_in,
        cex.lab = lab.font.size, cex.axis = ax.font.size, notch = TRUE)
grid(NULL,NULL, lty = 6)
boxplot(Energy~Method, data = all_metrics, 
        ylab = 'Energy', ylim = c(0, 1), col = colors_in, 
        medcol = colors_border, outcol = colors_in,
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, notch = TRUE)
abline(h = original.result["Energy"], col = ref_line_col,
       lty=2, lwd=3)
dev.off()

# Analysis of Age

pdf(file = "figures/Pembro_Age_Means.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.Age~Method, data = all_metrics, 
        ylab = 'Difference Means Age', ylim = c(0, max(all_metrics$M.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, notch = TRUE)
grid(NULL,NULL, lty = 6)
boxplot(M.Age~Method, data = all_metrics, 
        ylab = 'Difference Means Age', ylim = c(0, max(all_metrics$M.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, 
        col = colors_in,medcol = colors_border, outcol = colors_in, notch = TRUE)
abline(h = original.result["M.Age"], col = ref_line_col,
       lty=2, lwd=3)
dev.off()

pdf(file = "figures/Pembro_Age_SDs.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(SD.Age~Method, data = all_metrics, 
        ylab = 'Difference SD Age', ylim = c(0, max(all_metrics$SD.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, notch = TRUE)
grid(NULL,NULL, lty = 6)
boxplot(SD.Age~Method, data = all_metrics, 
        ylab = 'Difference SD Age', ylim = c(0, max(all_metrics$SD.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, 
        col = colors_in, medcol = colors_border, outcol = colors_in, notch = TRUE)
abline(h = original.result["SD.Age"], col = ref_line_col,
       lty=2, lwd=3)
dev.off()


# Analysis of PD.L1.MPS

pdf(file = "figures/Pembro_PD_Means.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.PD.L1.MPS~Method, data = all_metrics, 
        ylab = 'Difference Means PD-L1', ylim = c(0, max(all_metrics$M.PD.L1.MPS)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, notch = TRUE)
grid(NULL,NULL, lty = 6)
boxplot(M.PD.L1.MPS~Method, data = all_metrics, 
        ylab = 'Difference Means PD-L1', ylim = c(0, max(all_metrics$M.PD.L1.MPS)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, 
        col = colors_in, medcol = colors_border, outcol = colors_in, notch = TRUE)
abline(h = original.result["M.PD.L1.MPS"], col = ref_line_col,
       lty=2, lwd=3)
dev.off()

pdf(file = "figures/Pembro_PD_SDs.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(SD.PD.L1.MPS~Method, data = all_metrics, 
        ylab = 'Difference SD PD-L1', ylim = c(0, max(all_metrics$SD.PD.L1.MPS)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, notch = TRUE)
grid(NULL,NULL, lty = 6)
boxplot(SD.PD.L1.MPS~Method, data = all_metrics, 
        ylab = 'Difference SD PD-L1', ylim = c(0, max(all_metrics$SD.PD.L1.MPS)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, 
        col = colors_in, medcol = colors_border, outcol = colors_in, notch = TRUE)
abline(h = original.result["SD.PD.L1.MPS"], col = ref_line_col,
       lty=2, lwd=3)
dev.off()

# Analysis of Hemoglobin

pdf(file = "figures/Pembro_Hemoglobin_Means.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.Hemoglobin~Method, data = all_metrics, 
        ylab = 'Difference Means Hemoglobin', ylim = c(0, max(all_metrics$M.Hemoglobin)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, notch = TRUE)
grid(NULL,NULL, lty = 6)
boxplot(M.Hemoglobin~Method, data = all_metrics, 
        ylab = 'Difference Means Hemoglobin', ylim = c(0, max(all_metrics$M.Hemoglobin)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, 
        col = colors_in, medcol = colors_border, outcol = colors_in, notch = TRUE)
abline(h = original.result["M.Hemoglobin"], col = ref_line_col,
       lty=2, lwd=3)
dev.off()

pdf(file = "figures/Pembro_Hemoglobin_SDs.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(SD.Hemoglobin~Method, data = all_metrics, 
        ylab = 'Difference SD Hemoglobin', ylim = c(0, max(all_metrics$SD.Hemoglobin)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, notch = TRUE)
grid(NULL,NULL, lty = 6)
boxplot(SD.Hemoglobin~Method, data = all_metrics, 
        ylab = 'Difference SD Hemoglobin', ylim = c(0, max(all_metrics$SD.Hemoglobin)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, col = colors_in, 
        medcol = colors_border, outcol = colors_in, notch = TRUE)
abline(h = original.result["SD.Hemoglobin"], col = ref_line_col,
       lty=2, lwd=3)
dev.off()


pdf(file = "figures/Pembro_CG_Mean.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.CG~Method, data = all_metrics, yaxt = "n",
        ylab = 'Mean Correct Guess Probability', ylim = c(0,1),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, notch = TRUE)
axis(2, at = seq(0, 1, by = 0.25), cex.lab = lab.font.size, cex.axis = ax.font.size)
grid(NULL,NA, lty = 6, equilogs=F)
abline(h = 0.25, col = "lightgray", lty = "dotdash", lwd = par("lwd"))
abline(h = 0.5, col = "lightgray", lty = "dotdash", lwd = par("lwd"))
abline(h = 0.75, col = "lightgray", lty = "dotdash", lwd = par("lwd"))
abline(h = 1, col = "lightgray", lty = "dotdash", lwd = par("lwd"))
boxplot(M.CG~Method, data = all_metrics, yaxt = "n",
        ylab = 'Mean Correct Guess Probability',
        cex.lab = lab.font.size, cex.axis = ax.font.size, add = TRUE, 
        col = colors_in, medcol = colors_border, outcol = colors_in, notch = TRUE)
dev.off()


# Radar plot.
set.seed(99)
aggregate.metrics <- aggregate(all_metrics[,-1], by = list(all_metrics$Method), 
                               FUN = mean)
data.p <- aggregate.metrics[, c(2, 9, 10)] # aggregate.metrics[,-c(1, 11)]
row.names(data.p) <- aggregate.metrics[,1]
colnames(data.p) <- c("Balance in Group Sizes", "Energy", "Mean Correct Guess Prob.")

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data.p <- rbind(rep(2, ncol(data.p)), rep(0, ncol(data.p)), data.p)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
pdf(file = "figures/Pembro_Radar_Mean.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
radarchart( data.p[-c(1,2),]  , axistype=0 , maxmin=F,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=1.2 
)
legend(x=1, y=1, legend = rownames(data.p[-c(1,2),]), bty = "n", pch=20 , 
       col=colors_in , text.col = "black", cex=1.2, pt.cex=3)
dev.off()



