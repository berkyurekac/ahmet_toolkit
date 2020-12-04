# this script was created by Ahmet Can Berkyurek on 19/11/2019, to analyze the profile of mature piRNAs in rpb9 mutant. 
library(DESeq2)
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


setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/RBP9/5dependent")

#####
#introduce piRNA genes coverage at base resolution from wt, mj261 and rescue experiments. 
#####

wt1 <- read.delim("isa_10_5dependentAligned.sortedByCoord.out_adjusted.bed", header = FALSE, stringsAsFactors = FALSE)
wt3 <- read.delim("isa_12_5dependentAligned.sortedByCoord.out_adjusted.bed", header = FALSE, stringsAsFactors = FALSE)

mj261_1 <- read.delim("isa_13_5dependentAligned.sortedByCoord.out_adjusted.bed", header = FALSE, stringsAsFactors = FALSE)
mj261_2 <- read.delim("isa_14_5dependentAligned.sortedByCoord.out_adjusted.bed", header = FALSE, stringsAsFactors = FALSE)
mj261_3 <- read.delim("isa_15_5dependentAligned.sortedByCoord.out_adjusted.bed", header = FALSE, stringsAsFactors = FALSE)

rescue_1 <- read.delim("isa_16_5dependentAligned.sortedByCoord.out_adjusted.bed", header = FALSE, stringsAsFactors = FALSE)
rescue_2 <- read.delim("isa_17_5dependentAligned.sortedByCoord.out_adjusted.bed", header = FALSE, stringsAsFactors = FALSE)
rescue_3 <- read.delim("isa_18_5dependentAligned.sortedByCoord.out_adjusted.bed", header = FALSE, stringsAsFactors = FALSE)

wt1f <- 2634941	
wt3f <- 1813824	

mj261_1f <- 9000499	
mj261_2f <- 8124378	
mj261_3f <- 7396277	

rescue_1f <- 10924014	
rescue_2f <- 7661417
rescue_3f <- 3782503

# coberage is normalizes to reads per million
wt1$V8 <- (wt1$V8/wt1f)*1e6
wt3$V8 <- (wt3$V8/wt3f)*1e6

mj261_1$V8 <- (mj261_1$V8/mj261_1f)*1e6
mj261_2$V8 <- (mj261_2$V8/mj261_2f)*1e6
mj261_3$V8 <- (mj261_3$V8/mj261_3f)*1e6

rescue_1$V8 <- (rescue_1$V8/rescue_1f)*1e6
rescue_2$V8 <- (rescue_2$V8/rescue_2f)*1e6
rescue_3$V8 <- (rescue_3$V8/rescue_3f)*1e6

# combine replicates and plot density over a metagene piRNA
wtcombined <- data.frame(wt1$V4, wt1$V7, (wt1$V8+wt3$V8)/2)

mj261combined <- data.frame(mj261_1$V4, mj261_1$V7, (mj261_1$V8+mj261_2$V8+mj261_3$V8)/3)

rescuecombined <- data.frame(rescue_1$V4, rescue_1$V7, (rescue_1$V8+rescue_2$V8+rescue_3$V8)/3)


p <- ggplot(wtcombined, aes(x=wtcombined$wt1.V7, y=log10(wtcombined$X.wt1.V8...wt3.V8..2+1))) + geom_point(aes(colour = "black"), size=3) + theme_classic() + scale_y_continuous(limits = c(-1,4))

p1 <- ggplot(mj261combined, aes(x=mj261combined$mj261_1.V7, y=log10(mj261combined$X.mj261_1.V8...mj261_2.V8...mj261_3.V8..3+1))) + geom_point(aes(colour = "black"), size=3) + theme_classic() + scale_y_continuous(limits = c(-1,4))

p2 <- ggplot(rescuecombined, aes(x=rescuecombined$rescue_1.V7, y=log10(rescuecombined$X.rescue_1.V8...rescue_2.V8...rescue_3.V8..3+1))) + geom_point(aes(colour = "black"), size=3) + theme_classic() + scale_y_continuous(limits = c(-1,4))

ggsave("wt_piRNAprofiles_metagene.pdf", p, height = 10, width = 10, useDingbats=FALSE)
ggsave("mj261_piRNAprofiles_metagene.pdf", p1, height = 10, width = 10, useDingbats=FALSE)
ggsave("rescue_piRNAprofiles_metagene.pdf", p2, height = 10, width = 10, useDingbats=FALSE)



#new plot - density
wtcombined2 = wtcombined
colnames(wtcombined2) = c("name","x","y")
mj261combined2 = mj261combined
colnames(mj261combined2 ) = c("name","x","y")
rescuecombined2 = rescuecombined
colnames(rescuecombined2) = c("name","x","y")
wtcombined2$sample = "wt"
mj261combined2$sample = "mj261"
rescuecombined2$sample = "rescue"
combDF = rbind.data.frame(wtcombined2,mj261combined2,rescuecombined2)
g1 <- ggplot() + geom_smooth(data = combDF, aes(x = x, y = log10(y+1), colour = sample))+
  theme_classic()

ggsave("piRNAprofiles_metagene_all_smooth_3.pdf", g1, height = 15, width = 15, useDingbats=FALSE)


g <- ggplot() + geom_smooth(data = combDF[combDF$y != 0,], aes(x = x, y = log10(y+1), colour = sample))+
  theme_classic()

g <- ggplot() + geom_smooth(data = combDF[combDF$y != 0,], aes(x = x, y = log10(y+1), colour = sample))+
  theme_classic()

ggsave("piRNAprofiles_metagene_all_smooth_4.pdf", g, height = 15, width = 15, useDingbats=FALSE)


#new plot - heatmap -use reverse melting
#repeat the same on representative piRNAs


casted <- dcast(wtcombined2, wtcombined2$x)

#####
#####
#####

pirna <- read.csv("rpb9(mj261)vswt_allpiRNAs.csv", header = TRUE, stringsAsFactors = FALSE)
head(pirna)

sig <- subset(pirna, pirna$wt > pirna$mj261)
write.csv(sig, "sig_list.csv")
sig2 <- read.delim("")

pirna2 <- read.delim("1220952s_Additional_data_file_3.txt", header = TRUE, stringsAsFactors = FALSE)
head(pirna2)


a <- read.csv("lisa_ribozero_allgenes_padj.csv", header = TRUE, stringsAsFactors = FALSE)

head(a)

b <- read.delim("piRNA_targets_genes.tsv", header = FALSE, stringsAsFactors = FALSE)

head(b)


c <- subset(a, a$X %in% b$V1)
d <- subset(c, c$threshold == "A")
e <- subset(d, d$sig == "FDR<0.05")

f <- read.delim("piRNAtargets_top500.txt", header = FALSE, stringsAsFactors = FALSE)



p <- ggplot(wt1, aes(x=V7, y=V8)) + geom_point(aes(colour = "Black"), size=2) + scale_x_continuous(limits = c(0,32)) + theme_classic() + geom_smooth(method = "loess")









