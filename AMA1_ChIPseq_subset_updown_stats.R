# this script was created by Ahmet Can Berkyurek on 17/09/2020, to apply statistics on class I, II and III egenes from the AMA-1 ChIP-seq!
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library("genefilter")
library(pheatmap)
library(SummarizedExperiment)
library(gdata)
library(dplyr)
library(GenomicAlignments)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)


setwd("/Users/acanberkyurek/Dropbox (ericmiskalab)/RPB-9_discussion/Final Figures/EMBO/revision/bioinformatics/figure4/statistics_ChIP-seq/polyA")

##### 
# polyA data 
#####

##### 
# polyA upregulated
#####

##### 
# inroduce all replicates for class I, create a box plot for the TSS site, 100 bps up and down d
#####


wt_rep1 <- read.csv("wt_classI_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("wt_classI_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("wt_classI_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("mj261_classI_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("mj261_classI_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("mj261_classI_rep1.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame$meanwt <- rowMeans(frame[,2:4])
frame$meanmut <- rowMeans(frame[,5:7])
head(frame)
write.csv(frame, "AMA1_classI_Zscores.csv")

t.test(frame$meanwt, frame$meanmut)

ks.test(frame$meanwt, frame$meanmut)

wilcox.test(frame$meanwt, frame$meanmut) # pvalue 0.06353


a <- melt(frame)

g <- ggplot(a, aes(x=a$variable, y=a$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("AMA1_classI_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)


g2 <- ggplot(a, aes(x=a$variable, y=a$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("AMA1_classI_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)


##### 
# inroduce all replicates for class II, create a box plot for the TSS site, 100 bps up and down
#####


wt_rep1 <- read.csv("wt_classII_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("wt_classII_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("wt_classII_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("mj261_classII_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("mj261_classII_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("mj261_classII_rep1.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame2 <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame2) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame2$meanwt <- rowMeans(frame2[,2:4])
frame2$meanmut <- rowMeans(frame2[,5:7])
head(frame2)
write.csv(frame, "AMA1_classII_Zscores.csv")


t.test(frame2$meanwt, frame2$meanmut)

ks.test(frame2$meanwt, frame2$meanmut)

wilcox.test(frame2$meanwt, frame2$meanmut) # pvalue 2.062e-07


b <- melt(frame2)

g <- ggplot(b, aes(x=b$variable, y=b$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("AMA1_classII_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(b, aes(x=b$variable, y=b$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("AMA1_classII_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)



##### 
# inroduce all replicates for class III, create a box plot for the TSS site, 100 bps up and down
#####


wt_rep1 <- read.csv("wt_classIII_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("wt_classIII_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("wt_classIII_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("mj261_classIII_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("mj261_classIII_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("mj261_classIII_rep1.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame3 <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
frame3[,5:7] <- (frame[,5:7]/1.000255)
colnames(frame3) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame3$meanwt <- rowMeans(frame3[,2:4])
frame3$meanmut <- rowMeans(frame3[,5:7])
head(frame3)
write.csv(frame, "AMA1_classIII_Zscores.csv")


t.test(frame3$meanwt, frame3$meanmut)

ks.test(frame3$meanwt, frame3$meanmut)

wilcox.test(frame3$meanwt, frame3$meanmut) # 7.26e-05 pvalue 0.7765

c <- melt(frame3)

g <- ggplot(c, aes(x=c$variable, y=c$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("AMA1_classIII_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(c, aes(x=c$variable, y=c$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("AMA1_classIII_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)




setwd("/Users/acanberkyurek/Dropbox (ericmiskalab)/RPB-9_discussion/Final Figures/EMBO/revision/bioinformatics/figure4/statistics_ChIP-seq/polyA_downregulated")

##### 
# polyA downregulated
#####

##### 
# introduce all replicates for class I, create a box plot for the TSS site, 100 bps up and down d
#####


wt_rep1 <- read.csv("polyA_downregulated_classI_wt_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("polyA_downregulated_classI_wt_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("polyA_downregulated_classI_wt_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("polyA_downregulated_classI_mj261_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("polyA_downregulated_classI_mj261_rep2.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("polyA_downregulated_classI_mj261_rep3.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame$meanwt <- rowMeans(frame[,2:4])
frame$meanmut <- rowMeans(frame[,5:7])
head(frame)
write.csv(frame, "polyA_downregulated_AMA1_classI_Zscores.csv")

t.test(frame$meanwt, frame$meanmut)

ks.test(frame$meanwt, frame$meanmut)

wilcox.test(frame$meanwt, frame$meanmut) # pvalue 0.000391

##### p-value = 1.915e-05

a <- melt(frame)

g <- ggplot(a, aes(x=a$variable, y=a$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("polyA_downregulated_AMA1_classI_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)


g2 <- ggplot(a, aes(x=a$variable, y=a$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("polyA_ownregulated_AMA1_classI_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)


##### 
# inroduce all replicates for class II, create a box plot for the TSS site, 100 bps up and down
#####


wt_rep1 <- read.csv("polyA_downregulated_classII_wt_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("polyA_downregulated_classII_wt_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("polyA_downregulated_classII_wt_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("polyA_downregulated_classII_mj261_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("polyA_downregulated_classII_mj261_rep2.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("polyA_downregulated_classII_mj261_rep3.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame2 <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame2) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame2$meanwt <- rowMeans(frame2[,2:4])
frame2$meanmut <- rowMeans(frame2[,5:7])
head(frame2)
write.csv(frame2, "polyA_downregulated_AMA1_classII_Zscores.csv")


t.test(frame2$meanwt, frame2$meanmut)

ks.test(frame2$meanwt, frame2$meanmut)
##### p-value = 0.01042

wilcox.test(frame2$meanwt, frame2$meanmut) # pvalue 0.01356


b <- melt(frame2)

g <- ggplot(b, aes(x=b$variable, y=b$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("polyA_downregulated_AMA1_classII_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(b, aes(x=b$variable, y=b$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("polyA_downregulated_AMA1_classII_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)



##### 
# inroduce all replicates for class III, create a box plot for the TSS site, 100 bps up and down
#####



wt_rep1 <- read.csv("polyA_downregulated_classIII_wt_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("polyA_downregulated_classIII_wt_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("polyA_downregulated_classIII_wt_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("polyA_downregulated_classIII_mj261_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("polyA_downregulated_classIII_mj261_rep2.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("polyA_downregulated_classIII_mj261_rep3.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame3 <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
frame3[,5:7] <- (frame[,5:7]/1.000255)
colnames(frame3) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame3$meanwt <- rowMeans(frame3[,2:4])
frame3$meanmut <- rowMeans(frame3[,5:7])
head(frame3)
write.csv(frame3, "polyA_downregulated_AMA1_classIII_Zscores.csv")


t.test(frame3$meanwt, frame3$meanmut)

ks.test(frame3$meanwt, frame3$meanmut)
##### p-value = 0.567

wilcox.test(frame3$meanwt, frame3$meanmut) # pvalue 0.03089


c <- melt(frame3)

g <- ggplot(c, aes(x=c$variable, y=c$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("polyA_downregulated_AMA1_classIII_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(c, aes(x=c$variable, y=c$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("polyA_downregulated_AMA1_classIII_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)






setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Final Figures/EMBO/revision/bioinformatics/figure4/statistics_ChIP-seq/ribozero")

##### 
# ribozero data 
#####

##### 
# inroduce all replicates for class I, create a box plot for the TSS site, 100 bps up and down
#####


wt_rep1 <- read.csv("ribozero_wt_classI_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("ribozero_wt_classI_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("ribozero_wt_classI_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("ribozero_mj261_classI_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("ribozero_mj261_classI_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("ribozero_mj261_classI_rep1.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame$meanwt <- rowMeans(frame[,2:4])
frame$meanmut <- rowMeans(frame[,5:7])
head(frame)
write.csv(frame, "ribozero_AMA1_classI_Zscores.csv")

t.test(frame$meanwt, frame$meanmut)

ks.test(frame$meanwt, frame$meanmut)

##### p-value = 0.002342

a <- melt(frame)

g <- ggplot(a, aes(x=a$variable, y=a$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("ribozero_AMA1_classI_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)


g2 <- ggplot(a, aes(x=a$variable, y=a$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("ribozero_AMA1_classI_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)


##### 
# inroduce all replicates for class II, create a box plot for the TSS site, 100 bps up and down
#####


wt_rep1 <- read.csv("ribozero_wt_classII_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("ribozero_wt_classII_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("ribozero_wt_classII_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("ribozero_mj261_classII_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("ribozero_mj261_classII_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("ribozero_mj261_classII_rep1.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame2 <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame2) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame2$meanwt <- rowMeans(frame2[,2:4])
frame2$meanmut <- rowMeans(frame2[,5:7])
head(frame2)
write.csv(frame, "ribozero_AMA1_classII_Zscores.csv")


t.test(frame2$meanwt, frame2$meanmut)

ks.test(frame2$meanwt, frame2$meanmut)

#####  p-value = 0.0217

b <- melt(frame2)

g <- ggplot(b, aes(x=b$variable, y=b$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("ribozero_AMA1_classII_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(b, aes(x=b$variable, y=b$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("ribozero_AMA1_classII_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)



##### 
# inroduce all replicates for class III, create a box plot for the TSS site, 100 bps up and down
#####


wt_rep1 <- read.csv("ribozero_wt_classIII_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("ribozero_wt_classIII_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("ribozero_wt_classIII_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("ribozero_mj261_classIII_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("ribozero_mj261_classIII_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("ribozero_mj261_classIII_rep1.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame3 <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame3) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame3$meanwt <- rowMeans(frame3[,2:4])
frame3$meanmut <- rowMeans(frame3[,5:7])
frame3[,6:7] <- frame3[,6:7]/1.2
head(frame3)
write.csv(frame, "ribozero_AMA1_classIII_Zscores.csv")


t.test(frame3$meanwt, frame3$meanmut)

ks.test(frame3$meanwt, frame3$meanmut)

#####  p-value = 0.7912

c <- melt(frame3)

g <- ggplot(c, aes(x=c$variable, y=c$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("ribozero_AMA1_classIII_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(c, aes(x=c$variable, y=c$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("ribozero_AMA1_classIII_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)




##### 
# ribozero downregulated
#####

##### 
# introduce all replicates for class I, create a box plot for the TSS site, 100 bps up and down d
#####


wt_rep1 <- read.csv("polyA_downregulated_classI_wt_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("polyA_downregulated_classI_wt_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("polyA_downregulated_classI_wt_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("polyA_downregulated_classI_mj261_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("polyA_downregulated_classI_mj261_rep2.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("polyA_downregulated_classI_mj261_rep3.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame$meanwt <- rowMeans(frame[,2:4])
frame$meanmut <- rowMeans(frame[,5:7])
head(frame)
write.csv(frame, "polyA_downregulated_AMA1_classI_Zscores.csv")

t.test(frame$meanwt, frame$meanmut)

ks.test(frame$meanwt, frame$meanmut)
##### p-value = 1.915e-05

a <- melt(frame)

g <- ggplot(a, aes(x=a$variable, y=a$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("polyA_downregulated_AMA1_classI_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)


g2 <- ggplot(a, aes(x=a$variable, y=a$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("polyA_ownregulated_AMA1_classI_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)


##### 
# inroduce all replicates for class II, create a box plot for the TSS site, 100 bps up and down
#####


wt_rep1 <- read.csv("polyA_downregulated_classII_wt_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("polyA_downregulated_classII_wt_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("polyA_downregulated_classII_wt_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("polyA_downregulated_classII_mj261_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("polyA_downregulated_classII_mj261_rep2.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("polyA_downregulated_classII_mj261_rep3.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame2 <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
colnames(frame2) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame2$meanwt <- rowMeans(frame2[,2:4])
frame2$meanmut <- rowMeans(frame2[,5:7])
head(frame2)
write.csv(frame2, "polyA_downregulated_AMA1_classII_Zscores.csv")


t.test(frame2$meanwt, frame2$meanmut)

ks.test(frame2$meanwt, frame2$meanmut)
##### p-value = 0.01042

b <- melt(frame2)

g <- ggplot(b, aes(x=b$variable, y=b$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("polyA_downregulated_AMA1_classII_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(b, aes(x=b$variable, y=b$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("polyA_downregulated_AMA1_classII_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)



##### 
# inroduce all replicates for class III, create a box plot for the TSS site, 100 bps up and down
#####



wt_rep1 <- read.csv("polyA_downregulated_classIII_wt_rep1.csv", header = T, stringsAsFactors = F)

wt_rep2 <- read.csv("polyA_downregulated_classIII_wt_rep2.csv", header = T, stringsAsFactors = F)

wt_rep3 <- read.csv("polyA_downregulated_classIII_wt_rep3.csv", header = T, stringsAsFactors = F)

mj261_rep1 <- read.csv("polyA_downregulated_classIII_mj261_rep1.csv", header = T, stringsAsFactors = F)

mj261_rep2 <- read.csv("polyA_downregulated_classIII_mj261_rep2.csv", header = T, stringsAsFactors = F)

mj261_rep3 <- read.csv("polyA_downregulated_classIII_mj261_rep3.csv", header = T, stringsAsFactors = F)

head(wt_rep1)

frame3 <- data.frame(wt_rep1[,11:12], wt_rep2[,12], wt_rep3[,12], mj261_rep1[,12], mj261_rep2[,12], mj261_rep3[,12])
frame3[,5:7] <- (frame[,5:7]/1.000255)
colnames(frame3) <- c("order", "wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3")
frame3$meanwt <- rowMeans(frame3[,2:4])
frame3$meanmut <- rowMeans(frame3[,5:7])
head(frame3)
write.csv(frame3, "polyA_downregulated_AMA1_classIII_Zscores.csv")


t.test(frame3$meanwt, frame3$meanmut)

ks.test(frame3$meanwt, frame3$meanmut)
##### p-value = 0.000567

c <- melt(frame3)

g <- ggplot(c, aes(x=c$variable, y=c$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("wt_rep1", "wt_rep2", "wt_rep3", "mut_rep1", "mut_rep2", "mut_rep3"))

ggsave("polyA_downregulated_AMA1_classIII_allreplicates.pdf", g, height = 10, width = 10, useDingbats=FALSE)

g2 <- ggplot(c, aes(x=c$variable, y=c$value)) + geom_boxplot(outlier.colour = 'black') + theme_classic() + scale_x_discrete(limits=c("meanwt", "meanmut"))

ggsave("polyA_downregulated_AMA1_classIII_mean.pdf", g2, height = 10, width = 10, useDingbats=FALSE)


