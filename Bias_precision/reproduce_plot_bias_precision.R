###########################################################################
# Description:  This R script allows reproducing the plots in Sections 2.1
#               and 2.2 of Supplemental_simulations.pdf
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

library(grid)
library(ggplot2)
library(ggpubr)
# mainpath refers to the location of the folder 'Reproduce_simulations_results'
mainpath <- "C:/Users/Gabrielle/Google Drive/McGill - PhD/Thesis/dw-SM/JASA/Supplementary Material/Reproduce_simulations_results/"

# true blip parameters
psi1 <- c(0.1, 0.1)
psi2 <- c(-0.9, 0.6)


#### Visualize bias in scenarios 1-4, n=500, independent censoring 30% ####
secpath <- "Bias_precision/Independent30/"
setwd(paste(mainpath, secpath, sep = ""))
filelist <- list.files(pattern = ".*.txt")
datalist <-  lapply(filelist, FUN = read.table, header = TRUE)
datalist <- lapply(datalist, FUN = function(x){x[which(complete.cases(x)),]})

datalist1 <- rbind(datalist[[3]], datalist[[1]], datalist[[2]], datalist[[6]], datalist[[4]], datalist[[5]], datalist[[9]], datalist[[7]], datalist[[8]], datalist[[12]], datalist[[10]], datalist[[11]])
datalist2 <- rbind(datalist[[15]], datalist[[13]], datalist[[14]], datalist[[18]], datalist[[16]], datalist[[17]], datalist[[21]], datalist[[19]], datalist[[20]], datalist[[24]], datalist[[22]], datalist[[23]])
datalist3 <- rbind(datalist[[27]], datalist[[25]], datalist[[26]], datalist[[30]], datalist[[28]], datalist[[29]], datalist[[33]], datalist[[31]], datalist[[32]], datalist[[36]], datalist[[34]], datalist[[35]])
datalist4 <- rbind(datalist[[39]], datalist[[37]], datalist[[38]], datalist[[42]], datalist[[40]], datalist[[41]], datalist[[45]], datalist[[43]], datalist[[44]], datalist[[48]], datalist[[46]], datalist[[47]])

colnames(datalist1) <- c("psi0", "psi1")
np1 <- nrow(datalist1)/4
np1n1 <- nrow(datalist[[3]])
np1n2 <- nrow(datalist[[1]])
np1n3 <- nrow(datalist[[2]])
datalist1$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist1$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist1$scenario <- rep(1, nrow(datalist1))
datalist1$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist2) <- c("psi0", "psi1")
np1 <- nrow(datalist2)/4
np1n1 <- nrow(datalist[[15]])
np1n2 <- nrow(datalist[[13]])
np1n3 <- nrow(datalist[[14]])
datalist2$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist2$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist2$scenario <- rep(2, nrow(datalist2))
datalist2$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist3) <- c("psi0", "psi1")
np1 <- nrow(datalist3)/4
np1n1 <- nrow(datalist[[27]])
np1n2 <- nrow(datalist[[25]])
np1n3 <- nrow(datalist[[26]])
datalist3$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist3$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist3$scenario <- rep(3, nrow(datalist3))
datalist3$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist4) <- c("psi0", "psi1")
np1 <- nrow(datalist4)/4
np1n1 <- nrow(datalist[[39]])
np1n2 <- nrow(datalist[[37]])
np1n3 <- nrow(datalist[[38]])
datalist4$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist4$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist4$scenario <- rep(4, nrow(datalist4))
datalist4$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

alldata <- as.data.frame(rbind(datalist1, datalist2, datalist3, datalist4))

alldata1 <- alldata[alldata$stage == "1",]
alldata2 <- alldata[alldata$stage == "2",]
colnames(alldata1)[c(1,2)] <- c("psi10", "psi11")
colnames(alldata2)[c(1,2)] <- c("psi20", "psi21")
alldata1$scenario <- as.factor(alldata1$scenario)
alldata2$scenario <- as.factor(alldata2$scenario)
alldata1$n <- as.factor(alldata1$n)
alldata2$n <- as.factor(alldata2$n)
alldata1$method <- as.factor(alldata1$method)
alldata2$method <- as.factor(alldata2$method)

png(paste(mainpath, "Bias_precision/bias_precision_indep30_n500.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=10000, independent censoring 30% ####

png(paste(mainpath, "Bias_precision/bias_precision_indep30_n10000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()



#### Visualize bias in scenarios 1-4, n=500, dependent censoring 30% ####

secpath <- "Bias_precision/X1dependent30/"
setwd(paste(mainpath, secpath, sep = ""))
filelist <- list.files(pattern = ".*.txt")
datalist <-  lapply(filelist, FUN = read.table, header = TRUE)
datalist <- lapply(datalist, FUN = function(x){x[which(complete.cases(x)),]})

datalist1 <- rbind(datalist[[3]], datalist[[1]], datalist[[2]], datalist[[6]], datalist[[4]], datalist[[5]], datalist[[9]], datalist[[7]], datalist[[8]], datalist[[12]], datalist[[10]], datalist[[11]])
datalist2 <- rbind(datalist[[15]], datalist[[13]], datalist[[14]], datalist[[18]], datalist[[16]], datalist[[17]], datalist[[21]], datalist[[19]], datalist[[20]], datalist[[24]], datalist[[22]], datalist[[23]])
datalist3 <- rbind(datalist[[27]], datalist[[25]], datalist[[26]], datalist[[30]], datalist[[28]], datalist[[29]], datalist[[33]], datalist[[31]], datalist[[32]], datalist[[36]], datalist[[34]], datalist[[35]])
datalist4 <- rbind(datalist[[39]], datalist[[37]], datalist[[38]], datalist[[42]], datalist[[40]], datalist[[41]], datalist[[45]], datalist[[43]], datalist[[44]], datalist[[48]], datalist[[46]], datalist[[47]])

colnames(datalist1) <- c("psi0", "psi1")
np1 <- nrow(datalist1)/4
np1n1 <- nrow(datalist[[3]])
np1n2 <- nrow(datalist[[1]])
np1n3 <- nrow(datalist[[2]])
datalist1$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist1$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist1$scenario <- rep(1, nrow(datalist1))
datalist1$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist2) <- c("psi0", "psi1")
np1 <- nrow(datalist2)/4
np1n1 <- nrow(datalist[[15]])
np1n2 <- nrow(datalist[[13]])
np1n3 <- nrow(datalist[[14]])
datalist2$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist2$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist2$scenario <- rep(2, nrow(datalist2))
datalist2$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist3) <- c("psi0", "psi1")
np1 <- nrow(datalist3)/4
np1n1 <- nrow(datalist[[27]])
np1n2 <- nrow(datalist[[25]])
np1n3 <- nrow(datalist[[26]])
datalist3$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist3$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist3$scenario <- rep(3, nrow(datalist3))
datalist3$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist4) <- c("psi0", "psi1")
np1 <- nrow(datalist4)/4
np1n1 <- nrow(datalist[[39]])
np1n2 <- nrow(datalist[[37]])
np1n3 <- nrow(datalist[[38]])
datalist4$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist4$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist4$scenario <- rep(4, nrow(datalist4))
datalist4$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

alldata <- as.data.frame(rbind(datalist1, datalist2, datalist3, datalist4))

alldata1 <- alldata[alldata$stage == "1",]
alldata2 <- alldata[alldata$stage == "2",]
colnames(alldata1)[c(1,2)] <- c("psi10", "psi11")
colnames(alldata2)[c(1,2)] <- c("psi20", "psi21")
alldata1$scenario <- as.factor(alldata1$scenario)
alldata2$scenario <- as.factor(alldata2$scenario)
alldata1$n <- as.factor(alldata1$n)
alldata2$n <- as.factor(alldata2$n)
alldata1$method <- as.factor(alldata1$method)
alldata2$method <- as.factor(alldata2$method)

png(paste(mainpath, "Bias_precision/bias_precision_dep30_n500.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=1000, dependent censoring 30% ####

png(paste(mainpath, "Bias_precision/bias_precision_dep30_n1000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=10000, dependent censoring 30% ####

png(paste(mainpath, "Bias_precision/bias_precision_dep30_n10000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=500, independent censoring 60% ####

secpath <- "Bias_precision/Independent60/"
setwd(paste(mainpath, secpath, sep = ""))
filelist <- list.files(pattern = ".*.txt")
datalist <-  lapply(filelist, FUN = read.table, header = TRUE)
datalist <- lapply(datalist, FUN = function(x){x[which(complete.cases(x)),]})

datalist1 <- rbind(datalist[[3]], datalist[[1]], datalist[[2]], datalist[[6]], datalist[[4]], datalist[[5]], datalist[[9]], datalist[[7]], datalist[[8]], datalist[[12]], datalist[[10]], datalist[[11]])
datalist2 <- rbind(datalist[[15]], datalist[[13]], datalist[[14]], datalist[[18]], datalist[[16]], datalist[[17]], datalist[[21]], datalist[[19]], datalist[[20]], datalist[[24]], datalist[[22]], datalist[[23]])
datalist3 <- rbind(datalist[[27]], datalist[[25]], datalist[[26]], datalist[[30]], datalist[[28]], datalist[[29]], datalist[[33]], datalist[[31]], datalist[[32]], datalist[[36]], datalist[[34]], datalist[[35]])
datalist4 <- rbind(datalist[[39]], datalist[[37]], datalist[[38]], datalist[[42]], datalist[[40]], datalist[[41]], datalist[[45]], datalist[[43]], datalist[[44]], datalist[[48]], datalist[[46]], datalist[[47]])

colnames(datalist1) <- c("psi0", "psi1")
np1 <- nrow(datalist1)/4
np1n1 <- nrow(datalist[[3]])
np1n2 <- nrow(datalist[[1]])
np1n3 <- nrow(datalist[[2]])
datalist1$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist1$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist1$scenario <- rep(1, nrow(datalist1))
datalist1$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist2) <- c("psi0", "psi1")
np1 <- nrow(datalist2)/4
np1n1 <- nrow(datalist[[15]])
np1n2 <- nrow(datalist[[13]])
np1n3 <- nrow(datalist[[14]])
datalist2$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist2$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist2$scenario <- rep(2, nrow(datalist2))
datalist2$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist3) <- c("psi0", "psi1")
np1 <- nrow(datalist3)/4
np1n1 <- nrow(datalist[[27]])
np1n2 <- nrow(datalist[[25]])
np1n3 <- nrow(datalist[[26]])
datalist3$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist3$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist3$scenario <- rep(3, nrow(datalist3))
datalist3$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist4) <- c("psi0", "psi1")
np1 <- nrow(datalist4)/4
np1n1 <- nrow(datalist[[39]])
np1n2 <- nrow(datalist[[37]])
np1n3 <- nrow(datalist[[38]])
datalist4$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist4$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist4$scenario <- rep(4, nrow(datalist4))
datalist4$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

alldata <- as.data.frame(rbind(datalist1, datalist2, datalist3, datalist4))

alldata1 <- alldata[alldata$stage == "1",]
alldata2 <- alldata[alldata$stage == "2",]
colnames(alldata1)[c(1,2)] <- c("psi10", "psi11")
colnames(alldata2)[c(1,2)] <- c("psi20", "psi21")
alldata1$scenario <- as.factor(alldata1$scenario)
alldata2$scenario <- as.factor(alldata2$scenario)
alldata1$n <- as.factor(alldata1$n)
alldata2$n <- as.factor(alldata2$n)
alldata1$method <- as.factor(alldata1$method)
alldata2$method <- as.factor(alldata2$method)

png(paste(mainpath, "Bias_precision/bias_precision_indep60_n500.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=1000, independent censoring 60% ####

png(paste(mainpath, "Bias_precision/bias_precision_indep60_n1000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=10000, independent censoring 60% ####

png(paste(mainpath, "Bias_precision/bias_precision_indep60_n10000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=500, dependent censoring 60% ####

secpath <- "Bias_precision/X1dependent60/"
setwd(paste(mainpath, secpath, sep = ""))
filelist <- list.files(pattern = ".*.txt")
datalist <-  lapply(filelist, FUN = read.table, header = TRUE)
datalist <- lapply(datalist, FUN = function(x){x[which(complete.cases(x)),]})

datalist1 <- rbind(datalist[[3]], datalist[[1]], datalist[[2]], datalist[[6]], datalist[[4]], datalist[[5]], datalist[[9]], datalist[[7]], datalist[[8]], datalist[[12]], datalist[[10]], datalist[[11]])
datalist2 <- rbind(datalist[[15]], datalist[[13]], datalist[[14]], datalist[[18]], datalist[[16]], datalist[[17]], datalist[[21]], datalist[[19]], datalist[[20]], datalist[[24]], datalist[[22]], datalist[[23]])
datalist3 <- rbind(datalist[[27]], datalist[[25]], datalist[[26]], datalist[[30]], datalist[[28]], datalist[[29]], datalist[[33]], datalist[[31]], datalist[[32]], datalist[[36]], datalist[[34]], datalist[[35]])
datalist4 <- rbind(datalist[[39]], datalist[[37]], datalist[[38]], datalist[[42]], datalist[[40]], datalist[[41]], datalist[[45]], datalist[[43]], datalist[[44]], datalist[[48]], datalist[[46]], datalist[[47]])

colnames(datalist1) <- c("psi0", "psi1")
np1 <- nrow(datalist1)/4
np1n1 <- nrow(datalist[[3]])
np1n2 <- nrow(datalist[[1]])
np1n3 <- nrow(datalist[[2]])
datalist1$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist1$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist1$scenario <- rep(1, nrow(datalist1))
datalist1$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist2) <- c("psi0", "psi1")
np1 <- nrow(datalist2)/4
np1n1 <- nrow(datalist[[15]])
np1n2 <- nrow(datalist[[13]])
np1n3 <- nrow(datalist[[14]])
datalist2$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist2$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist2$scenario <- rep(2, nrow(datalist2))
datalist2$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist3) <- c("psi0", "psi1")
np1 <- nrow(datalist3)/4
np1n1 <- nrow(datalist[[27]])
np1n2 <- nrow(datalist[[25]])
np1n3 <- nrow(datalist[[26]])
datalist3$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist3$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist3$scenario <- rep(3, nrow(datalist3))
datalist3$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

colnames(datalist4) <- c("psi0", "psi1")
np1 <- nrow(datalist4)/4
np1n1 <- nrow(datalist[[39]])
np1n2 <- nrow(datalist[[37]])
np1n3 <- nrow(datalist[[38]])
datalist4$stage <- as.factor(rep(c(rep(1, np1), rep(2, np1)), 2))
datalist4$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 4)
datalist4$scenario <- rep(4, nrow(datalist4))
datalist4$method <- c(rep("dWSurv", np1*2), rep("HNW", np1*2))

alldata <- as.data.frame(rbind(datalist1, datalist2, datalist3, datalist4))

alldata1 <- alldata[alldata$stage == "1",]
alldata2 <- alldata[alldata$stage == "2",]
colnames(alldata1)[c(1,2)] <- c("psi10", "psi11")
colnames(alldata2)[c(1,2)] <- c("psi20", "psi21")
alldata1$scenario <- as.factor(alldata1$scenario)
alldata2$scenario <- as.factor(alldata2$scenario)
alldata1$n <- as.factor(alldata1$n)
alldata2$n <- as.factor(alldata2$n)
alldata1$method <- as.factor(alldata1$method)
alldata2$method <- as.factor(alldata2$method)

png(paste(mainpath, "Bias_precision/bias_precision_dep60_n500.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=1000, dependent censoring 60% ####

png(paste(mainpath, "Bias_precision/bias_precision_dep60_n1000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=10000, dependent censoring 60% ####

png(paste(mainpath, "Bias_precision/bias_precision_dep60_n10000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi10, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi11, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi20, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi21, fill = method)) + geom_boxplot() + theme_minimal() + scale_fill_manual(values = c("grey35", "grey69")) + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=500, time-varying dependent censoring 30% ####

secpath <- "Bias_precision/X1X2dependent30/"
setwd(paste(mainpath, secpath, sep = ""))
filelist <- list.files(pattern = ".*.txt")
datalist <-  lapply(filelist, FUN = read.table, header = TRUE)
datalist <- lapply(datalist, FUN = function(x){x[which(complete.cases(x)),]})

datalist1 <- rbind(datalist[[3]], datalist[[1]], datalist[[2]], datalist[[6]], datalist[[4]], datalist[[5]])
datalist2 <- rbind(datalist[[9]], datalist[[7]], datalist[[8]], datalist[[12]], datalist[[10]], datalist[[11]])
datalist3 <- rbind(datalist[[15]], datalist[[13]], datalist[[14]], datalist[[18]], datalist[[16]], datalist[[17]])
datalist4 <- rbind(datalist[[21]], datalist[[19]], datalist[[20]], datalist[[24]], datalist[[22]], datalist[[23]])

colnames(datalist1) <- c("psi0", "psi1")
np1 <- nrow(datalist1)/2
np1n1 <- nrow(datalist[[3]])
np1n2 <- nrow(datalist[[1]])
np1n3 <- nrow(datalist[[2]])
datalist1$stage <- as.factor(c(rep(1, np1), rep(2, np1)))
datalist1$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 2)
datalist1$scenario <- rep(1, nrow(datalist1))

colnames(datalist2) <- c("psi0", "psi1")
np1 <- nrow(datalist2)/2
np1n1 <- nrow(datalist[[9]])
np1n2 <- nrow(datalist[[7]])
np1n3 <- nrow(datalist[[8]])
datalist2$stage <- as.factor(c(rep(1, np1), rep(2, np1)))
datalist2$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 2)
datalist2$scenario <- rep(2, nrow(datalist2))

colnames(datalist3) <- c("psi0", "psi1")
np1 <- nrow(datalist3)/2
np1n1 <- nrow(datalist[[15]])
np1n2 <- nrow(datalist[[13]])
np1n3 <- nrow(datalist[[14]])
datalist3$stage <- as.factor(c(rep(1, np1), rep(2, np1)))
datalist3$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 2)
datalist3$scenario <- rep(3, nrow(datalist3))

colnames(datalist4) <- c("psi0", "psi1")
np1 <- nrow(datalist4)/2
np1n1 <- nrow(datalist[[21]])
np1n2 <- nrow(datalist[[19]])
np1n3 <- nrow(datalist[[20]])
datalist4$stage <- as.factor(c(rep(1, np1), rep(2, np1)))
datalist4$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 2)
datalist4$scenario <- rep(4, nrow(datalist4))

alldata <- as.data.frame(rbind(datalist1, datalist2, datalist3, datalist4))

alldata1 <- alldata[alldata$stage == "1",]
alldata2 <- alldata[alldata$stage == "2",]
colnames(alldata1)[c(1,2)] <- c("psi10", "psi11")
colnames(alldata2)[c(1,2)] <- c("psi20", "psi21")
alldata1$scenario <- as.factor(alldata1$scenario)
alldata2$scenario <- as.factor(alldata2$scenario)
alldata1$n <- as.factor(alldata1$n)
alldata2$n <- as.factor(alldata2$n)

png(paste(mainpath, "Bias_precision/bias_precision_X1X2dep30_n500.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi10)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi11)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi20)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi21)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=1000, time-varying dependent censoring 30% ####

png(paste(mainpath, "Bias_precision/bias_precision_X1X2dep30_n1000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi10)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi11)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi20)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi21)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=1000, time-varying dependent censoring 30% ####

png(paste(mainpath, "Bias_precision/bias_precision_X1X2dep30_n10000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi10)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi11)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi20)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi21)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=500, time-varying dependent censoring 60% ####

secpath <- "Bias_precision/X1X2dependent60/"
setwd(paste(mainpath, secpath, sep = ""))
filelist <- list.files(pattern = ".*.txt")
datalist <-  lapply(filelist, FUN = read.table, header = TRUE)
datalist <- lapply(datalist, FUN = function(x){x[which(complete.cases(x)),]})

datalist1 <- rbind(datalist[[3]], datalist[[1]], datalist[[2]], datalist[[6]], datalist[[4]], datalist[[5]])
datalist2 <- rbind(datalist[[9]], datalist[[7]], datalist[[8]], datalist[[12]], datalist[[10]], datalist[[11]])
datalist3 <- rbind(datalist[[15]], datalist[[13]], datalist[[14]], datalist[[18]], datalist[[16]], datalist[[17]])
datalist4 <- rbind(datalist[[21]], datalist[[19]], datalist[[20]], datalist[[24]], datalist[[22]], datalist[[23]])

colnames(datalist1) <- c("psi0", "psi1")
np1 <- nrow(datalist1)/2
np1n1 <- nrow(datalist[[3]])
np1n2 <- nrow(datalist[[1]])
np1n3 <- nrow(datalist[[2]])
datalist1$stage <- as.factor(c(rep(1, np1), rep(2, np1)))
datalist1$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 2)
datalist1$scenario <- rep(1, nrow(datalist1))

colnames(datalist2) <- c("psi0", "psi1")
np1 <- nrow(datalist2)/2
np1n1 <- nrow(datalist[[9]])
np1n2 <- nrow(datalist[[7]])
np1n3 <- nrow(datalist[[8]])
datalist2$stage <- as.factor(c(rep(1, np1), rep(2, np1)))
datalist2$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 2)
datalist2$scenario <- rep(2, nrow(datalist2))

colnames(datalist3) <- c("psi0", "psi1")
np1 <- nrow(datalist3)/2
np1n1 <- nrow(datalist[[15]])
np1n2 <- nrow(datalist[[13]])
np1n3 <- nrow(datalist[[14]])
datalist3$stage <- as.factor(c(rep(1, np1), rep(2, np1)))
datalist3$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 2)
datalist3$scenario <- rep(3, nrow(datalist3))

colnames(datalist4) <- c("psi0", "psi1")
np1 <- nrow(datalist4)/2
np1n1 <- nrow(datalist[[21]])
np1n2 <- nrow(datalist[[19]])
np1n3 <- nrow(datalist[[20]])
datalist4$stage <- as.factor(c(rep(1, np1), rep(2, np1)))
datalist4$n <- rep(c(rep(500, np1n1), rep(1000, np1n2), rep(10000, np1n3)), 2)
datalist4$scenario <- rep(4, nrow(datalist4))

alldata <- as.data.frame(rbind(datalist1, datalist2, datalist3, datalist4))

alldata1 <- alldata[alldata$stage == "1",]
alldata2 <- alldata[alldata$stage == "2",]
colnames(alldata1)[c(1,2)] <- c("psi10", "psi11")
colnames(alldata2)[c(1,2)] <- c("psi20", "psi21")
alldata1$scenario <- as.factor(alldata1$scenario)
alldata2$scenario <- as.factor(alldata2$scenario)
alldata1$n <- as.factor(alldata1$n)
alldata2$n <- as.factor(alldata2$n)

png(paste(mainpath, "Bias_precision/bias_precision_X1X2dep60_n500.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi10)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "500",], aes(x = scenario, y = psi11)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi20)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "500",], aes(x = scenario, y = psi21)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=1000, time-varying dependent censoring 60% ####

png(paste(mainpath, "Bias_precision/bias_precision_X1X2dep60_n1000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi10)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "1000",], aes(x = scenario, y = psi11)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi20)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "1000",], aes(x = scenario, y = psi21)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

#### Visualize bias in scenarios 1-4, n=1000, time-varying dependent censoring 60% ####

png(paste(mainpath, "Bias_precision/bias_precision_X1X2dep60_n10000.png", sep = ""), height = 650, width = 800)
p10 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi10)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[10])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-0.3, 0.9), breaks=c(-0.3, 0, 0.3, 0.6, 0.9))
p11 <- ggplot(alldata1[alldata1$n == "10000",], aes(x = scenario, y = psi11)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi1[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab(expression(hat(psi)[11])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(-1.1, 0.7), breaks=c(-1, -0.6, -0.2, 0.2, 0.6))
p20 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi20)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[1], size = 0.8, linetype="dashed") + theme(legend.title = element_blank()) + ylab(expression(hat(psi)[20])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75)) + scale_y_continuous(limits = c(-1.9, -0.2), breaks = c(-1.9, -1.5, -1.1, -0.7, -0.3))
p21 <- ggplot(alldata2[alldata2$n == "10000",], aes(x = scenario, y = psi21)) + geom_boxplot() + theme_minimal() + geom_hline(yintercept = psi2[2], size = 0.8, linetype="dashed") + theme(legend.title = element_blank(), axis.title.y = element_text(margin = margin(t = 0, r = 8.5, b = 0, l = 0))) + ylab(expression(hat(psi)[21])) + stat_summary(fun.y = mean, geom = "point", shape = 8, size = 3, position = position_dodge(width = 0.75))  + scale_y_continuous(limits = c(0.1, 1.3), breaks = c(0.1, 0.4, 0.7, 1, 1.3))
ggarrange(p10, p11, p20, p21, ncol = 2, nrow = 2, heights = c(2, 2.2), labels = c("A)", "B)", "C)", "D)"), common.legend = TRUE)
dev.off()

