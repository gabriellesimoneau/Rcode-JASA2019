###########################################################################
# Description:  This R script allows reproducing the plots in Section 4 of
#               Supplemental_simulations.pdf
#
# User instructions: The folder 'Reproduce_simulations_results' should be
#               saved on the user's computer. The variable mainpath should
#               be replaced by the location of that folder on the user's
#               computer.
#
# Important note: The following script depends on how the simulation 
#               results in .txt are sorted in their corresponding folders. 
#               The .txt files should be sorted by name.
#
# Last updated: July 26, 2018
###########################################################################

library(ggplot2)
# mainpath: location of the folder "Reproduce_simulations_results"
mainpath <- ""
secpath <- "Comparison_survival_time/Results/"
setwd(paste(mainpath, secpath, sep = ""))


#### Visualize survival time distribution, n=500 and 1000, independent censoring 30% ####
d500 <- read.table("indep30_500.txt")
d1000 <- read.table("indep30_1000.txt")

dall <- as.data.frame(c(as.matrix(d500), as.matrix(d1000)))
colnames(dall) <- "survtime"
dall$method <- as.factor(rep(rep(c("True optimal", "dWSurv", "HNW", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"), each = 10000), 2))
dall$n <- as.factor(rep(c(500, 1000), each = 70000))
dall$method <- factor(dall$method, levels = c("True optimal", "dWSurv", "HNW", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"))
dall <- dall[which(dall$method %in% c("dWSurv", "HNW") | dall$n == "500"),]

png(paste(mainpath, "Comparison_survival_time/compare_indep30.png", sep = ""), width = 643, height = 375)
ggplot(dall, aes(x= method, y = survtime, fill = n)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("white", "grey69")) + xlab("Method") + ylab("Survival time")
dev.off()

#### Visualize survival time distribution, n=500 and 1000, independent censoring 60% ####
d500 <- read.table("indep60_500.txt")
d1000 <- read.table("indep60_1000.txt")

dall <- as.data.frame(c(as.matrix(d500), as.matrix(d1000)))
colnames(dall) <- "survtime"
dall$method <- as.factor(rep(rep(c("True optimal", "dWSurv", "HNW", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"), each = 10000), 2))
dall$n <- as.factor(rep(c(500, 1000), each = 70000))
dall$method <- factor(dall$method, levels = c("True optimal", "dWSurv", "HNW", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"))
dall <- dall[which(dall$method %in% c("dWSurv", "HNW") | dall$n == "500"),]

png(paste(mainpath, "Comparison_survival_time/compare_indep60.png", sep = ""), width = 643, height = 375)
ggplot(dall, aes(x= method, y = survtime, fill = n)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("white", "grey69")) + xlab("Method") + ylab("Survival time")
dev.off()

#### Visualize survival time distribution, n=500 and 1000, baseline dependent censoring 30% ####
d500 <- read.table("X1dep30_500.txt")
d1000 <- read.table("X1dep30_1000.txt")

dall <- as.data.frame(c(as.matrix(d500), as.matrix(d1000)))
colnames(dall) <- "survtime"
dall$method <- as.factor(rep(rep(c("True optimal", "dWSurv", "HNW", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"), each = 10000), 2))
dall$n <- as.factor(rep(c(500, 1000), each = 70000))
dall$method <- factor(dall$method, levels = c("True optimal", "dWSurv", "HNW", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"))
dall <- dall[which(dall$method %in% c("dWSurv", "HNW") | dall$n == "500"),]

png(paste(mainpath, "Comparison_survival_time/compare_X1dep30.png", sep = ""), width = 643, height = 375)
ggplot(dall, aes(x= method, y = survtime, fill = n)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("white", "grey69")) + xlab("Method") + ylab("Survival time")
dev.off()

#### Visualize survival time distribution, n=500 and 1000, baseline dependent censoring 60% ####
d500 <- read.table("X1dep60_500.txt")
d1000 <- read.table("X1dep60_1000.txt")

dall <- as.data.frame(c(as.matrix(d500), as.matrix(d1000)))
colnames(dall) <- "survtime"
dall$method <- as.factor(rep(rep(c("True optimal", "dWSurv", "HNW", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"), each = 10000), 2))
dall$n <- as.factor(rep(c(500, 1000), each = 70000))
dall$method <- factor(dall$method, levels = c("True optimal", "dWSurv", "HNW", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"))
dall <- dall[which(dall$method %in% c("dWSurv", "HNW") | dall$n == "500"),]

png(paste(mainpath, "Comparison_survival_time/compare_X1dep60.png", sep = ""), width = 643, height = 375)
ggplot(dall, aes(x= method, y = survtime, fill = n)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("white", "grey69")) + xlab("Method") + ylab("Survival time")
dev.off()

#### Visualize survival time distribution, n=500 and 1000, time-varying dependent censoring 30% ####
d500 <- read.table("X1X2dep30_500.txt")
d1000 <- read.table("X1X2dep30_1000.txt")

dall <- as.data.frame(c(as.matrix(d500), as.matrix(d1000)))
colnames(dall) <- "survtime"
dall$method <- as.factor(rep(rep(c("True optimal", "dWSurv", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"), each = 10000), 2))
dall$n <- as.factor(rep(c(500, 1000), each = 60000))
dall$method <- factor(dall$method, levels = c("True optimal", "dWSurv", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"))
dall <- dall[which(dall$method == "dWSurv" | dall$n == "500"),]

png(paste(mainpath, "Comparison_survival_time/compare_X1X2dep30.png", sep = ""), width = 643, height = 375)
ggplot(dall, aes(x= method, y = survtime, fill = n)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("white", "grey69")) + xlab("Method") + ylab("Survival time")
dev.off()

#### Visualize survival time distribution, n=500 and 1000, time-varying dependent censoring 60% ####
d500 <- read.table("X1X2dep60_500.txt")
d1000 <- read.table("X1X2dep60_1000.txt")

dall <- as.data.frame(c(as.matrix(d500), as.matrix(d1000)))
colnames(dall) <- "survtime"
dall$method <- as.factor(rep(rep(c("True optimal", "dWSurv", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"), each = 10000), 2))
dall$n <- as.factor(rep(c(500, 1000), each = 60000))
dall$method <- factor(dall$method, levels = c("True optimal", "dWSurv", "A1=A2=0", "A1=0,A2=1", "A1=1,A2=0", "A1=A2=1"))
dall <- dall[which(dall$method == "dWSurv" | dall$n == "500"),]

png(paste(mainpath, "Comparison_survival_time/compare_X1X2dep60.png", sep = ""), width = 643, height = 375)
ggplot(dall, aes(x= method, y = survtime, fill = n)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("white", "grey69")) + xlab("Method") + ylab("Survival time")
dev.off()