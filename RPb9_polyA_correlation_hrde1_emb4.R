# this script was created by Ahmet Can Berkyurek on 06/04/2020, to run correlation between RPB-9 polyA RNA-seq and hrde-1, emb-4 total RNA-seq!
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


#####
# introduce rpb-9 mut diff. genes, then introduce hrde-1 and emb-4 diff. genes
#####

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_correlation/")

rpb9 <- read.csv("mj261_ribozero_diff_osc.csv", header = TRUE, stringsAsFactors = FALSE)
head(rpb9)



rpb9sig <- subset(rpb9, !(rpb9$threshold) == "C")

hrde1 <- read.csv("hrde1_allgenes_osc.csv", header = TRUE, stringsAsFactors = FALSE)
head(hrde1)

emb4 <- read.csv("emb4_diff.csv", header = TRUE, stringsAsFactors = FALSE)
head(emb4)

hrde1 <- read.csv("prg-1(n4357)_diff_osc_mongomery_gonaddissected.csv", header = TRUE, stringsAsFactors = FALSE)
head(prg1)

hrde1 <- read.csv("prg-1(alper)_diff_osc.csv", header = TRUE, stringsAsFactors = FALSE)
head(hrde1)

germline <- read.csv("germline_expressed.csv", header = TRUE, stringsAsFactors = FALSE)

soma <- read.delim("soma_genes_reinke.bed", header = F)

rpb9soma <- subset(rpb9, rpb9$X %in% soma$V1)

hrde1soma <- subset(hrde1, hrde1$X %in% soma$V1)

#####
# merge data sets - rpb9 and prg-1
#####

merged <- merge(rpb9soma, hrde1soma, by.x = "X", by.y = "X")
head(merged)
merged <- data.frame(merged)

merged <- subset(merged, merged$X %in% germline$genes)


p <- ggplot(merged, aes(x=merged$log2FoldChange.x, y=merged$log2FoldChange.y))  + geom_point(color = "Grey", size=3)  + scale_x_continuous(limits = c(-5,5)) + scale_y_continuous(limits = c(-5,5)) + geom_vline(xintercept = 0, colour = "grey") + geom_abline(slope =0.36848, intercept = 0.016178  ) + theme_classic() + xlab("rpb-9(mj261)_ribo") + ylab("prg-1_alper_ribo")

ggsave("rpb9ribo_prg1_alper_correlation_somatic.pdf", p, height = 10, width = 10, useDingbats=FALSE)



reg = lm(merged$log2FoldChange.x ~merged$log2FoldChange.y)
summary(reg)  ##### 0.03405   0.008002   31.33   <2e-16 ***


rpb9up <- subset(rpb9, rpb9$threshold == "A")
rpb9upmerged <- subset(merged, merged$X %in% rpb9up$X)
rpb9down <- subset(rpb9, rpb9$threshold == "B")
rpb9downmerged <- subset(merged, merged$X %in% rpb9down$X)

hrde1up <- subset(hrde1, hrde1$threshold == "A")
hrde1down <- subset(hrde1, hrde1$threshold == "B")


rpb9uphrde1 <- subset(rpb9up, rpb9up$X %in% hrde1up$X)
write.csv(rpb9uphrde1, "mj261ribo_up_csr1cat_common.csv")
rpb9uphrde1merged <- subset(merged, merged$X %in% rpb9uphrde1$X)
rpb9downhrde1 <- subset(rpb9down, rpb9down$X %in% hrde1down$X)
write.csv(rpb9downhrde1, "mj261ribo_down_csr1cat_common.csv")
rpb9downhrde1merged <- subset(merged, merged$X %in% rpb9downhrde1$X)

sperm <- read.delim("sperm_genes_reinke.bed", header = FALSE)
sperm_merged <- subset(merged, merged$X %in% sperm$V1)
soma <- read.delim("soma_genes_reinke.bed", header = FALSE)
soma_merged <- subset(merged, merged$X %in% soma$V1)

sensor <- subset(merged, merged$X == "piRNASensor")

somamerged <- subset(rpb9uphrde1merged, rpb9uphrde1merged$X %in% soma$V1)


p_rpb9_prg1 <- ggplot(merged, aes(x=merged$log2FoldChange.x, y=merged$log2FoldChange.y))  + geom_point(color = "Grey", size=3) + geom_point(data = somamerged, aes(x=somamerged$log2FoldChange.x, y=somamerged$log2FoldChange.y), color="yellow", size=7) + geom_point(data=sensor, aes(x=sensor$log2FoldChange.x, y=sensor$log2FoldChange.y), color="green", size=7.5) + geom_point(data=rpb9uphrde1merged, aes(x=rpb9uphrde1merged$log2FoldChange.x, y=rpb9uphrde1merged$log2FoldChange.y), color="red", size=4.5) + geom_point(data=rpb9downhrde1merged, aes(x=rpb9downhrde1merged$log2FoldChange.x, y=rpb9downhrde1merged$log2FoldChange.y), color="blue", size=4.5) + geom_point(data=rpb9upmerged, aes(x=rpb9upmerged$log2FoldChange.x, y=rpb9upmerged$log2FoldChange.y), color="black", size=1.5) + geom_point(data=rpb9downmerged, aes(x=rpb9downmerged$log2FoldChange.x, y=rpb9downmerged$log2FoldChange.y), color="black", size=1.5) + scale_x_continuous(limits = c(-12,12)) + scale_y_continuous(limits = c(-12,12)) + geom_vline(xintercept = 0, colour = "grey") + geom_abline(slope = 0.0399, intercept = 0.016178  ) + theme_classic() + xlab("rpb-9(mj261)_ribo") + ylab("prg-1_alper_ribo")

ggsave("rpb9ribo_csr1_catmut_correlation_gonad.pdf", p_rpb9_prg1, height = 10, width = 10, useDingbats=FALSE)


p_rpb9_prg1 <- ggplot(merged, aes(x=merged$log2FoldChange.x, y=merged$log2FoldChange.y))  + geom_point(color = "Grey", size=3) + geom_point(data=rpb9uphrde1merged, aes(x=rpb9uphrde1merged$log2FoldChange.x, y=rpb9uphrde1merged$log2FoldChange.y), color="red", size=4.5) + geom_point(data=rpb9downhrde1merged, aes(x=rpb9downhrde1merged$log2FoldChange.x, y=rpb9downhrde1merged$log2FoldChange.y), color="blue", size=4.5) + geom_point(data=rpb9upmerged, aes(x=rpb9upmerged$log2FoldChange.x, y=rpb9upmerged$log2FoldChange.y), color="black", size=1.5) + geom_point(data=rpb9downmerged, aes(x=rpb9downmerged$log2FoldChange.x, y=rpb9downmerged$log2FoldChange.y), color="black", size=1.5) + scale_x_continuous(limits = c(-12,12)) + scale_y_continuous(limits = c(-12,12)) + geom_vline(xintercept = 0, colour = "grey") + 
  geom_abline(slope =0.334074, intercept = 0.008002 ) + theme_classic() + xlab("rpb-9(mj261)_polyA") + ylab("prg-1_gonad")

ggsave("./all_replicates/rpb9polyA_prg1(n4357)_correlation_gonad_allreplicates.pdf", p_rpb9_prg1, height = 10, width = 10, useDingbats=FALSE)



#####
# merge data sets - rpb9 and hrde1
#####

merged <- merge(rpb9, hrde1, by.x = "X", by.y = "X")
head(merged)
merged <- data.frame(merged)


# merged <- subset(merged, merged$X %in% germline$genes)


reg = lm(merged$log2FoldChange.x ~merged$log2FoldChange.y)
summary(reg)  ##### R2 0.09812, R2  0.03576 


rpb9up <- subset(rpb9, rpb9$threshold == "A")
rpb9upmerged <- subset(merged, merged$X %in% rpb9up$X)
rpb9down <- subset(rpb9, rpb9$threshold == "B")
rpb9downmerged <- subset(merged, merged$X %in% rpb9down$X)

hrde1up <- subset(hrde1, hrde1$threshold == "A")
hrde1down <- subset(hrde1, hrde1$threshold == "B")


rpb9uphrde1 <- subset(rpb9up, rpb9up$X %in% hrde1up$X)
rpb9uphrde1merged <- subset(merged, merged$X %in% rpb9uphrde1$X)
piRNASensor <- subset(rpb9uphrde1merged, rpb9uphrde1merged$X == "piRNASensor")
rpb9downhrde1 <- subset(rpb9down, rpb9down$X %in% hrde1down$X)
rpb9downhrde1merged <- subset(merged, merged$X %in% rpb9downhrde1$X)

sperm <- read.delim("sperm_genes_reinke.bed", header = FALSE)
sperm_merged <- subset(merged, merged$X %in% sperm$V1)
soma <- read.delim("soma_genes_reinke.bed", header = FALSE)
soma_merged <- subset(merged, merged$X %in% soma$V1)


p_rpb9_hrde1 <- ggplot(merged, aes(x=merged$log2FoldChange.x, y=merged$log2FoldChange.y))  + geom_point(color = "Grey", size=3)  + geom_point(data=piRNASensor, aes(x=piRNASensor$log2FoldChange.x, y=piRNASensor$log2FoldChange.y), color="Green", size=7.5) + geom_point(data=rpb9uphrde1merged, aes(x=rpb9uphrde1merged$log2FoldChange.x, y=rpb9uphrde1merged$log2FoldChange.y), color="red", size=4.5) + geom_point(data=rpb9downhrde1merged, aes(x=rpb9downhrde1merged$log2FoldChange.x, y=rpb9downhrde1merged$log2FoldChange.y), color="blue", size=4.5) + geom_point(data=rpb9upmerged, aes(x=rpb9upmerged$log2FoldChange.x, y=rpb9upmerged$log2FoldChange.y), color="black", size=1.5) + geom_point(data=rpb9downmerged, aes(x=rpb9downmerged$log2FoldChange.x, y=rpb9downmerged$log2FoldChange.y), color="black", size=1.5) + scale_x_continuous(limits = c(-12,12)) + scale_y_continuous(limits = c(-12,12)) + geom_vline(xintercept = 0, colour = "grey") + 
  geom_abline(slope =0.1812, intercept = 0.01277 ) + theme_classic() + xlab("rpb-9(mj261)_polyA_allreplicates") + ylab("hrde1_ribo")

ggsave("rpb9polyA_hrde1_correlation_allreplicates.pdf", p_rpb9_hrde1, height = 10, width = 10, useDingbats=FALSE)


#####
# merge data sets - rpb9 and hrde1 - shuffle the data frame to see if the overlap is not random.
#####

# set a random seeed so that it is reproducible. 
set.seed(1001)

merged2 <- merged

merged2$log2FoldChange.x <- sample(merged2$log2FoldChange.x)
merged2$log2FoldChange.y <- sample(merged2$log2FoldChange.y)


reg = lm(merged2$log2FoldChange.x ~merged2$log2FoldChange.y)
summary(reg)  ##### R2 -8.574e-07, R2 7.299e-05 ,  -4.039e-05 


rpb9up <- subset(rpb9, rpb9$threshold == "A")
rpb9upmerged <- subset(merged2, merged2$X %in% rpb9up$X)
rpb9down <- subset(rpb9, rpb9$threshold == "B")
rpb9downmerged <- subset(merged2, merged2$X %in% rpb9down$X)

hrde1up <- subset(hrde1, hrde1$threshold == "A")
hrde1down <- subset(hrde1, hrde1$threshold == "B")


rpb9uphrde1 <- subset(rpb9up, rpb9up$X %in% hrde1up$X)
rpb9uphrde1merged <- subset(merged2, merged2$X %in% rpb9uphrde1$X)
piRNASensor <- subset(rpb9uphrde1merged, rpb9uphrde1merged$X == "piRNASensor")
rpb9downhrde1 <- subset(rpb9down, rpb9down$X %in% hrde1down$X)
rpb9downhrde1merged <- subset(merged2, merged2$X %in% rpb9downhrde1$X)


p_rpb9_hrde1 <- ggplot(merged2, aes(x=merged2$log2FoldChange.x, y=merged2$log2FoldChange.y))  + geom_point(color = "Grey", size=3)  + geom_point(data=piRNASensor, aes(x=piRNASensor$log2FoldChange.x, y=piRNASensor$log2FoldChange.y), color="Green", size=7.5) + geom_point(data=rpb9uphrde1merged, aes(x=rpb9uphrde1merged$log2FoldChange.x, y=rpb9uphrde1merged$log2FoldChange.y), color="red", size=4.5)  + geom_point(data=rpb9upmerged, aes(x=rpb9upmerged$log2FoldChange.x, y=rpb9upmerged$log2FoldChange.y), color="black", size=1.5) + geom_point(data=rpb9downmerged, aes(x=rpb9downmerged$log2FoldChange.x, y=rpb9downmerged$log2FoldChange.y), color="black", size=1.5) + scale_x_continuous(limits = c(-12,12)) + scale_y_continuous(limits = c(-12,12)) + geom_vline(xintercept = 0, colour = "grey") + 
  geom_abline(slope =-0.013, intercept = 0.01277 ) + theme_classic() + xlab("rpb-9(mj261)_polyA_allreplicates") + ylab("hrde1_ribo")

ggsave("rpb9_prg1_correlation_shuffled.pdf", p_rpb9_hrde1, height = 10, width = 10, useDingbats=FALSE)


#####
# merge data sets - rpb9 and prg1 - shuffle the data frame to see if the overlap is not random.
#####

# set a random seeed so that it is reproducible. 
set.seed(1001)

merged2 <- merged

merged2$log2FoldChange.x <- sample(merged2$log2FoldChange.x)
merged2$log2FoldChange.y <- sample(merged2$log2FoldChange.y)


reg = lm(merged2$log2FoldChange.x ~merged2$log2FoldChange.y)
summary(reg)  ##### 0.05733   0.01277   31.33   <2e-16 ***


rpb9up <- subset(rpb9, rpb9$threshold == "A")
rpb9upmerged <- subset(merged2, merged2$X %in% rpb9up$X)
rpb9down <- subset(rpb9, rpb9$threshold == "B")
rpb9downmerged <- subset(merged2, merged2$X %in% rpb9down$X)

hrde1up <- subset(hrde1, hrde1$threshold == "A")
hrde1down <- subset(hrde1, hrde1$threshold == "B")


rpb9uphrde1 <- subset(rpb9up, rpb9up$X %in% hrde1up$X)
rpb9uphrde1merged <- subset(merged2, merged2$X %in% rpb9uphrde1$X)
piRNASensor <- subset(rpb9uphrde1merged, rpb9uphrde1merged$X == "piRNASensor")
rpb9downhrde1 <- subset(rpb9down, rpb9down$X %in% hrde1down$X)
rpb9downhrde1merged <- subset(merged2, merged2$X %in% rpb9downhrde1$X)


p_rpb9_hrde1 <- ggplot(merged2, aes(x=merged2$log2FoldChange.x, y=merged2$log2FoldChange.y))  + geom_point(color = "Grey", size=3)  + geom_point(data=piRNASensor, aes(x=piRNASensor$log2FoldChange.x, y=piRNASensor$log2FoldChange.y), color="Green", size=7.5) + geom_point(data=rpb9uphrde1merged, aes(x=rpb9uphrde1merged$log2FoldChange.x, y=rpb9uphrde1merged$log2FoldChange.y), color="red", size=4.5)  + geom_point(data=rpb9upmerged, aes(x=rpb9upmerged$log2FoldChange.x, y=rpb9upmerged$log2FoldChange.y), color="black", size=1.5) + geom_point(data=rpb9downmerged, aes(x=rpb9downmerged$log2FoldChange.x, y=rpb9downmerged$log2FoldChange.y), color="black", size=1.5) + scale_x_continuous(limits = c(-12,12)) + scale_y_continuous(limits = c(-12,12)) + geom_vline(xintercept = 0, colour = "grey") + 
  geom_abline(slope =0.02025, intercept = 0.01277 ) + theme_classic() + xlab("rpb-9(mj261)_polyA_allreplicates") + ylab("hrde1_ribo")

ggsave("rpb9_hrde1_correlation_shuffled.pdf", p_rpb9_hrde1, height = 10, width = 10, useDingbats=FALSE)



# + geom_point(data=rpb9downhrde1merged, aes(x=rpb9downhrde1merged$log2FoldChange.x, y=rpb9downhrde1merged$log2FoldChange.y), color="blue", size=4.5)

#####
# merge data sets - rpb9 and emb-4
#####

merged2 <- merge(rpb9, emb4, by.x = "X", by.y = "X")
head(merged2)
merged <- data.frame(merged2)


reg = lm(merged2$log2FoldChange.x ~merged2$log2FoldChange.y)
summary(reg)  #####0.1897   0.006389   62.97   <2e-16 ***


rpb9up <- subset(rpb9, rpb9$threshold == "A")
rpb9upmerged <- subset(merged2, merged2$X %in% rpb9up$X)
rpb9down <- subset(rpb9, rpb9$threshold == "B")
rpb9downmerged <- subset(merged2, merged2$X %in% rpb9down$X)

hrde1up <- subset(emb4, emb4$threshold == "A")
hrde1down <- subset(emb4, emb4$threshold == "B")


rpb9uphrde1 <- subset(rpb9up, rpb9up$X %in% hrde1up$X)
rpb9uphrde1merged <- subset(merged2, merged2$X %in% rpb9uphrde1$X)
piRNASensor <- subset(rpb9uphrde1merged, rpb9uphrde1merged$X == "piRNASensor")
rpb9downhrde1 <- subset(rpb9down, rpb9down$X %in% hrde1down$X)
write.csv(rpb9uphrde1, "mj271_emb4_common.csv")
rpb9downhrde1merged <- subset(merged2, merged2$X %in% rpb9downhrde1$X)



p_rpb9_emb4 <- ggplot(merged2, aes(x=merged2$log2FoldChange.x, y=merged2$log2FoldChange.y))  + geom_point(color = "Grey", size=3) + geom_point(data=piRNASensor, aes(x=piRNASensor$log2FoldChange.x, y=piRNASensor$log2FoldChange.y), color="Green", size=7.5) + geom_point(data=rpb9uphrde1merged, aes(x=rpb9uphrde1merged$log2FoldChange.x, y=rpb9uphrde1merged$log2FoldChange.y), color="red", size=4.5) + geom_point(data=rpb9downhrde1merged, aes(x=rpb9downhrde1merged$log2FoldChange.x, y=rpb9downhrde1merged$log2FoldChange.y), color="blue", size=4.5) + geom_point(data=rpb9upmerged, aes(x=rpb9upmerged$log2FoldChange.x, y=rpb9upmerged$log2FoldChange.y), color="black", size=1.5) + geom_point(data=rpb9downmerged, aes(x=rpb9downmerged$log2FoldChange.x, y=rpb9downmerged$log2FoldChange.y), color="black", size=1.5) + scale_x_continuous(limits = c(-12,12)) + scale_y_continuous(limits = c(-12,12)) + geom_vline(xintercept = 0, colour = "grey") + 
  geom_abline(slope =  0.1995, intercept = 0.006389 ) + theme_classic() + xlab("rpb-9(mj261)_polyA_allreplicates") + ylab("emb4_ribo")

ggsave("rpb9polyA_emb4_correlation.pdf", p_rpb9_emb4, height = 10, width = 10, useDingbats=FALSE)




