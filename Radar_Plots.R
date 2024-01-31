# Library
library(fmsb)
library(RColorBrewer)
library(scales)

# Create data: note in High school for several students
set.seed(99)
aggregate.metrics <- aggregate(all_metrics[,-1], by = list(all_metrics$Method), FUN= mean)

data.p <- aggregate.metrics[,-c(1, 11)]
colnames(data.p) <- c("Balance", "M.Age", "M.PD.L1", "M.Hemoglobin", "SD.Age", "SD.PD.L1", 
                      "SD.Hemoglobin", "Energy", "M.CG")
row.names(data.p) <- aggregate.metrics[,1]
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
#data.p <- rbind(apply(data.p, MARGIN = 2, FUN = max), apply(data.p, MARGIN = 2, FUN = max), data.p)
data.p <- rbind(rep(1.5, ncol(data.p)), rep(0, ncol(data.p)), data.p)

#radarchart( data.p )  

coul <- brewer.pal(5, "Set1")
colors_border <- coul
colors_in <- alpha(coul,0.3)

# If you remove the 2 first lines, the function compute the max and min of each variable with the available data:
pdf(file = "Radar_test.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
radarchart( data.p[-c(1,2),]  , axistype=0 , maxmin=F,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8, 
            #custom labels
            vlcex=1.2 
)
legend(x=1, y=1, legend = rownames(data.p[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.2, pt.cex=3)
dev.off()

pdf(file = "Pembro_Age_Means_Test.pdf", width = 11, height = 8) # defaults to 7 x 7 inches
boxplot(M.Age~Method, data = all_metrics, 
        ylab = 'Difference Means Age', ylim = c(0, max(all_metrics$M.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, outpch = 19)
grid(NULL,NULL, lty = 6)
boxplot(M.Age~Method, data = all_metrics, 
        ylab = 'Difference Means Age', ylim = c(0, max(all_metrics$M.Age)),
        cex.lab = lab.font.size, cex.axis = ax.font.size, col = colors_in, 
        medcol = colors_border, outcol = colors_in, outpch = 19,
        add = TRUE)
abline(h = original.result["M.Age"], col = 'red',
       lty=2, lwd=3)
dev.off()
