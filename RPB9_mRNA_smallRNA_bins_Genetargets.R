# this script was created by Ahmet Can Berkyurek on 10/04/2020, to correlate mRNA and smallRNa expression in bins and analayze gene targets
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
library(DESeq2)

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/others/smallRNA_bins")

#load the normalized data

data <- read.csv("RPB9polyA_mRNA_smallRNA_bins_log10TPM.csv", header = TRUE, stringsAsFactors = FALSE)

head(data)


# first, subset with piRNA targets and other arganoute targets


piRNAtop <- read.delim("./gene_targets/piRNAtargets.txt", header = FALSE)
piRNA2MM<- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs/2MM_piRNatargets.txt", header = TRUE)
CSR1 <- read.delim("./gene_targets/CSR-1.txt", header = FALSE)
HRDE1 <- read.delim("./gene_targets/HRDE1.txt", header = FALSE)
ERGO1 <- read.delim("./gene_targets/ERGO-1.txt", header = FALSE)
ALG34 <- read.delim("./gene_targets/ALG-34.txt", header = FALSE)
germlineRNAi <- read.delim("./gene_targets/germline.txt", header = FALSE)
somaRNAi <- read.delim("./gene_targets/soma.txt", header = FALSE)
wago <- read.delim("./gene_targets/WAGO.txt", header = FALSE)


data_wago <- subset(data, data$genes %in% wago$V1)
freqwago <- as.data.frame(table(data_wago$label))
freqwago$abundance <- (freqwago$Freq/sum(freqwago$Freq))*100
write.csv(data_wago, "data_wago.csv")

data_piRNA <- subset(data, data$genes %in% piRNAtop$V1)
freqpiRNA <- as.data.frame(table(data_piRNA$label))
freqpiRNA$abundance <- (freqpiRNA$Freq/sum(freqpiRNA$Freq))*100
write.csv(freqpiRNA, "./gene_targets/Freq_piRNA_polyA.csv")
write.csv(data_piRNA, "data_piRNA.csv")

data_piRNA2MM <- subset(data, data$genes %in% piRNA2MM$Gene.stable.ID)
freqpiRNA2MM <- as.data.frame(table(data_piRNA2MM$label))
freqpiRNA2MM$abundance <- (freqpiRNA2MM$Freq/sum(freqpiRNA2MM$Freq))*100
write.csv(data_piRNA2MM, "data_piRNA2MM.csv")

data_CSR1 <- subset(data, data$genes %in% CSR1$V1)
freqCSR1 <- as.data.frame(table(data_CSR1$label))
freqCSR1$abundance <- (freqCSR1$Freq/sum(freqCSR1$Freq))*100
write.csv(freqCSR1, "./gene_targets/Freq_CSR1_polyA.csv")
write.csv(data_CSR1, "data_CSR1.csv")

data_HRDE1 <- subset(data, data$genes %in% HRDE1$V1)
freqHRDE1 <- as.data.frame(table(data_HRDE1$label))
freqHRDE1$abundance <- (freqHRDE1$Freq/sum(freqHRDE1$Freq))*100
write.csv(freqHRDE1, "./gene_targets/Freq_HRDE1_polyA.csv")
write.csv(data_HRDE1, "data_HRDE1.csv")

data_ERGO1 <- subset(data, data$genes %in% ERGO1$V1)
freqERGO1 <- as.data.frame(table(data_ERGO1$label))
freqERGO1$abundance <- (freqERGO1$Freq/sum(freqERGO1$Freq))*100
write.csv(data_ERGO1, "data_ERGO1.csv")


data_ALG34 <- subset(data, data$genes %in% ALG34$V1)
freqALG34 <- as.data.frame(table(data_ALG34$label))
freqALG34$abundance <- (freqALG34$Freq/sum(freqALG34$Freq))*100
write.csv(data_ALG34, "data_ALG34.csv")


data_germlineRNAi <- subset(data, data$genes %in% germlineRNAi$V1)
freqgermlineRNAi <- as.data.frame(table(data_germlineRNAi$label))
freqgermlineRNAi$abundance <- (freqgermlineRNAi$Freq/sum(freqgermlineRNAi$Freq))*100
write.csv(data_germlineRNAi, "data_germlineRNAi.csv")


data_somaRNAi <- subset(data, data$genes %in% somaRNAi$V1)
freqsomaRNAi <- as.data.frame(table(data_somaRNAi$label))
freqsomaRNAi$abundance <- (freqsomaRNAi$Freq/sum(freqsomaRNAi$Freq))*100
write.csv(data_somaRNAi, "data_somaRNAi.csv")


seq <- as.vector(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U"))

piRNAx <- data.frame(c(0,0,0,0,0.2645503,0,0 , 0.2645503  ,0.7936508 , 0.5291005 , 1.8518519 , 2.3809524 , 3.1746032  ,3.1746032,  2.1164021 , 4.7619048,  4.4973545 , 6.3492063, 16.6666667 ,53.1746032),seq)

CSR1x <- data.frame(c(0.2590674,0,0,0,0,  0.5181347,0,0,0,0,0,  0.7772021 , 1.0362694,  4.1450777 , 8.8082902, 15.2849741 ,15.8031088, 22.2797927 ,21.7616580 , 9.3264249),seq)

HRDE1x <- data.frame(c(0,0,0,0,0,0,0,0,0.2197802,0,0,  0.2197802 , 0.6593407,  0.2197802 , 0.6593407,  0.6593407 , 0.8791209 , 1.9780220, 12.5274725, 81.9780220),seq)

ERGO1x <- data.frame(c(0,0,0,0,0,0,0,0,0,0,0, 2.439024,  2.439024  ,2.439024 , 2.439024,0,0,  2.439024 , 7.317073 ,80.487805),seq)

wagox <- data.frame(c(0,0,0,0,0,0,0,0,0,0,0.09250694,0,0.09250694,0.09250694 , 0.27752081 , 0.83256244 , 3.42275671 , 8.14061055, 19.98149861, 67.06753006),seq)

ALG34x <- data.frame(c(0,0.2493766,0,0,0,0,0,0,0.7481297, 0.7481297 , 1.4962594,  8.2294264 , 9.7256858, 13.9650873 ,17.4563591 ,14.2144638 , 4.9875312 ,11.2219451  ,9.7256858,  7.2319202),seq)


all <- data.frame(piRNAx$seq, piRNAx$c.0..0..0..0..0.2645503..0..0..0.2645503..0.7936508..0.5291005.., wagox$c.0..0..0..0..0.09624639..0..0.09624639..0..0.09624639..0..0.38498556..  , HRDE1x$c.0..0..0..0..0..0..0..0..0.2197802..0..0..0.2197802..0.6593407.., ERGO1x$c.0..0..0..0..0..0..0..0..0..0..0..2.439024..2.439024..2.439024.., ALG34x$c.0.2518892..0.2518892..0.2518892..0.2518892..3.5264484..4.0302267..,CSR1x$c.0.2590674..0..0..0..0..0.5181347..0..0..0..0..0..0.7772021..)
colnames(all) <- c("label", "piRNA", "WAGO", "HRDE-1", "ERGO-1","ALG-3/4","CSR-1")
print(all)

all2 <- melt(all)

p <- ggplot(all2, aes(x=label, y=value, fill=variable)) + geom_bar(position="dodge", stat = "identity") + theme_classic() + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins log10 (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets.pdf", p, height = 10, width = 15, useDingbats=FALSE)



all <- data.frame(piRNAx$seq, piRNAx$c.0..0..0..0..0.2645503..0..0..0.2645503..0.7936508..0.5291005.. , HRDE1x$c.0..0..0..0..0..0..0..0..0.2197802..0..0..0.2197802..0.6593407.., CSR1x$c.0.2590674..0..0..0..0..0.5181347..0..0..0..0..0..0.7772021..)
colnames(all) <- c("label", "piRNA", "HRDE-1","CSR-1")
print(all)

all2 <- melt(all)

p <- ggplot(all2, aes(x=label, y=value, fill=variable)) + geom_bar(position="dodge", stat = "identity") + theme_classic() + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins log10 (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets_onlyHRDE1CSR1.pdf", p, height = 10, width = 15, useDingbats=FALSE)


#load the normalized data

data <- read.csv("RPB9polyA_mRNA_smallRNA_bins_log10TPM_last4bins.csv", header = TRUE, stringsAsFactors = FALSE)

head(data)


# then, re-cluster the last 4 bins. this might give a more sensitive information. 


piRNAtop <- read.delim("./gene_targets/piRNAtargets.txt", header = FALSE)
CSR1 <- read.delim("./gene_targets/CSR-1.txt", header = FALSE)
ERGO1 <- read.delim("./gene_targets/ERGO-1.txt", header = FALSE)
ALG34 <- read.delim("./gene_targets/ALG-34.txt", header = FALSE)
germlineRNAi <- read.delim("./gene_targets/germline.txt", header = FALSE)
somaRNAi <- read.delim("./gene_targets/soma.txt", header = FALSE)
wago <- read.delim("./gene_targets/WAGO.txt", header = FALSE)


data_wago <- subset(data, data$genes %in% wago$V1)
freqwago <- as.data.frame(table(data_wago$label))
freqwago$abundance <- (freqwago$Freq/sum(freqwago$Freq))*100
write.csv(data_wago, "data_wagosub.csv")

data_piRNA <- subset(data, data$genes %in% piRNAtop$V1)
freqpiRNA <- as.data.frame(table(data_piRNA$label))
freqpiRNA$abundance <- (freqpiRNA$Freq/sum(freqpiRNA$Freq))*100
write.csv(data_piRNA, "data_piRNAsub.csv")

data_CSR1 <- subset(data, data$genes %in% CSR1$V1)
freqCSR1 <- as.data.frame(table(data_CSR1$label))
freqCSR1$abundance <- (freqCSR1$Freq/sum(freqCSR1$Freq))*100
write.csv(data_CSR1, "data_CSR1sub.csv")


data_ERGO1 <- subset(data, data$genes %in% ERGO1$V1)
freqERGO1 <- as.data.frame(table(data_ERGO1$label))
freqERGO1$abundance <- (freqERGO1$Freq/sum(freqERGO1$Freq))*100
write.csv(data_ERGO1, "data_ERGO1sub.csv")


data_ALG34 <- subset(data, data$genes %in% ALG34$V1)
freqALG34 <- as.data.frame(table(data_ALG34$label))
freqALG34$abundance <- (freqALG34$Freq/sum(freqALG34$Freq))*100
write.csv(data_ALG34, "data_ALG34sub.csv")


data_germlineRNAi <- subset(data, data$genes %in% germlineRNAi$V1)
freqgermlineRNAi <- as.data.frame(table(data_germlineRNAi$label))
freqgermlineRNAi$abundance <- (freqgermlineRNAi$Freq/sum(freqgermlineRNAi$Freq))*100
write.csv(data_germlineRNAi, "data_germlineRNAisub.csv")


data_somaRNAi <- subset(data, data$genes %in% somaRNAi$V1)
freqsomaRNAi <- as.data.frame(table(data_somaRNAi$label))
freqsomaRNAi$abundance <- (freqsomaRNAi$Freq/sum(freqsomaRNAi$Freq))*100
write.csv(data_somaRNAi, "data_somaRNAisub.csv")


seq <- as.vector(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U"))

piRNAx <- data.frame(c(0.7874016,  1.5748031 , 0.5249344 , 1.0498688 , 2.3622047,  1.3123360,  1.3123360 , 1.0498688  ,1.5748031,  2.8871391 , 2.0997375,  1.8372703 , 3.1496063,  3.4120735 , 5.2493438 , 6.0367454, 11.0236220, 15.2230971, 18.6351706, 18.8976378),seq)

CSR1x <- data.frame(c(0,0.2631579 , 0.7894737 , 2.3684211 , 1.8421053 , 3.4210526,  3.4210526 , 5.0000000 , 7.8947368,  7.8947368  ,7.6315789 , 6.3157895, 10.5263158,  8.6842105 ,11.0526316 , 7.6315789,  7.8947368  ,5.5263158 , 1.8421053,0),seq)

ERGO1x <- data.frame(c(0,0,1.960784,0,3.921569,0,0,0,0,0,1.960784,3.921569,0,1.960784,  1.960784 , 3.921569,  3.921569 ,17.647059, 19.607843, 39.215686),seq)

wagox <- data.frame(c(0,0,0,0, 0.09624639 ,0, 0.09624639 ,0, 0.09624639 ,0, 0.38498556,  0.48123195,  1.25120308  ,2.21366699  ,2.98363811,  5.00481232  ,5.38979788,  9.43214629 ,20.01924928 ,52.55052936),seq)

ALG34x <- data.frame(c(0.2518892 , 0.2518892,  0.2518892 , 0.2518892 , 3.5264484 , 4.0302267 , 6.2972292,  9.5717884 , 8.3123426 , 8.5642317 ,11.8387909 ,10.3274559,  5.5415617 , 2.2670025  ,6.5491184  ,4.0302267, 4.0302267 , 5.2896725  ,4.7858942 , 4.0302267),seq)


all <- data.frame(piRNAx$seq, piRNAx$c.0..0..0..0..0.2645503..0..0..0.2645503..0.7936508..0.5291005.., wagox$c.0..0..0..0..0.09624639..0..0.09624639..0..0.09624639..0..0.38498556..  , HRDE1x$c.0..0..0..0..0..0..0..0..0.2197802..0..0..0.2197802..0.6593407.., ERGO1x$c.0..0..0..0..0..0..0..0..0..0..0..2.439024..2.439024..2.439024.., ALG34x$c.0.2518892..0.2518892..0.2518892..0.2518892..3.5264484..4.0302267..,CSR1x$c.0.2590674..0..0..0..0..0.5181347..0..0..0..0..0..0.7772021..)
colnames(all) <- c("label", "piRNA", "WAGO", "HRDE-1", "ERGO-1","ALG-3/4","CSR-1")
print(all)

all2 <- melt(all)

p <- ggplot(all2, aes(x=label, y=value, fill=variable)) + geom_bar(position="dodge", stat = "identity") + theme_classic() + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins log10 (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets.pdf", p, height = 10, width = 15, useDingbats=FALSE)


# then, repeat this for upregulated and ChIP-seq CLassI, II, III genes

polyA <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA/all_replicates/mj261_polyA_diff_osc_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)

up <- polyA[polyA$threshold =="A", ]

up <- na.omit(up)

down <- polyA[polyA$threshold =="B", ]

down <- na.omit(down)


data <- read.csv("RPB9polyA_mRNA_smallRNA_bins_log10TPM.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)

data_cluster <- read.csv("RPB9polyA_mRNA_smallRNA_bins_log10TPM_last4bins.csv", header = TRUE, stringsAsFactors = FALSE)
head(data_cluster)


classI <- read.delim("upclassI.tsv", header = FALSE)
classII <- read.delim("upclassII.tsv", header = FALSE)
classIII <- read.delim("upclassIII.tsv", header = FALSE)

data_up <- subset(data, data$genes %in% up$genes)
data_down <- subset(data, data$genes %in% down$genes)
data_upcluster <- subset(data_cluster, data_cluster$genes %in% up$genes)
data_downcluster <- subset(data_cluster, data_cluster$genes %in% down$genes)
data_classI <- subset(data, data$genes %in% classI$V1)
data_classII <- subset(data, data$genes %in% classII$V1)
data_classIII <- subset(data, data$genes %in% classIII$V1)

data_classI2 <- subset(data_cluster, data_cluster$genes %in% classI$V1)
data_classII2 <- subset(data_cluster, data_cluster$genes %in% classII$V1)
data_classIII2 <- subset(data_cluster, data_cluster$genes %in% classIII$V1)


freqdown <- as.data.frame(table(data_down$label))
freqdown$abundance <- (freqdown$Freq/sum(freqdown$Freq))*100
write.csv(data_down, "data_down.csv")

p <- ggplot(freqdown, aes(x=Var1, y=abundance)) + geom_bar( stat = "identity", color = "Black",width =0.25) + theme_classic() + scale_x_discrete(limits = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U")) + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets_downregulatedgenes.pdf", p, height = 10, width = 15, useDingbats=FALSE)


freqdown2 <- as.data.frame(table(data_downcluster$label))
freqdown2$abundance <- (freqdown2$Freq/sum(freqdown2$Freq))*100
write.csv(data_downcluster, "data_down_cluster.csv")


p <- ggplot(freqdown2, aes(x=Var1, y=abundance)) + geom_bar( stat = "identity", color = "Black",width =0.25) + theme_classic() + scale_x_discrete(limits = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U")) + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets_downregulatedgenes_last4clusters.pdf", p, height = 10, width = 15, useDingbats=FALSE)



frequp <- as.data.frame(table(data_up$label))
frequp$abundance <- (frequp$Freq/sum(frequp$Freq))*100
write.csv(frequp, "./gene_targets/Freq_upregulated_polyA.csv")
write.csv(data_up, "data_up.csv")

p <- ggplot(frequp, aes(x=Var1, y=abundance)) + geom_bar( stat = "identity", color = "Black",width =0.25) + theme_classic() + scale_x_discrete(limits = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U")) + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets_upregulatedgenes.pdf", p, height = 10, width = 15, useDingbats=FALSE)


frequp2 <- as.data.frame(table(data_upcluster$label))
frequp2$abundance <- (frequp2$Freq/sum(frequp2$Freq))*100
write.csv(data_upcluster, "data_up_cluster.csv")


p <- ggplot(frequp2, aes(x=Var1, y=abundance)) + geom_bar( stat = "identity", color = "Black",width =0.25) + theme_classic() + scale_x_discrete(limits = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U")) + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets_upregulatedgenes_last4clusters.pdf", p, height = 10, width = 15, useDingbats=FALSE)


freqclassI <- as.data.frame(table(data_classI$label))
freqclassI$abundance <- (freqclassI$Freq/sum(freqclassI$Freq))*100
write.csv(freqclassI, "./gene_targets/Freq_ClassI_polyA.csv")
write.csv(data_classI, "data_classI.csv")

freqclassII <- as.data.frame(table(data_classII$label))
freqclassII$abundance <- (freqclassII$Freq/sum(freqclassII$Freq))*100
write.csv(freqclassII, "./gene_targets/Freq_ClassII_polyA.csv")
write.csv(data_classII, "data_classII.csv")


freqclassIII <- as.data.frame(table(data_classIII$label))
freqclassIII$abundance <- (freqclassIII$Freq/sum(freqclassIII$Freq))*100
write.csv(freqclassIII, "./gene_targets/Freq_ClassIII_polyA.csv")
write.csv(data_classIII, "data_classIII.csv")

freqclassI2 <- as.data.frame(table(data_classI2$label))
freqclassI2$abundance <- (freqclassI2$Freq/sum(freqclassI2$Freq))*100
write.csv(data_classI2, "data_classIcluster.csv")

freqclassII2 <- as.data.frame(table(data_classII2$label))
freqclassII2$abundance <- (freqclassII2$Freq/sum(freqclassII2$Freq))*100
write.csv(data_classII2, "data_classIIcluster.csv")


freqclassIII2 <- as.data.frame(table(data_classIII2$label))
freqclassIII2$abundance <- (freqclassIII2$Freq/sum(freqclassIII2$Freq))*100
write.csv(data_classIII2, "data_classIIIcluster.csv")


seq <- as.vector(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U"))

upx <- data.frame(c(0.7874016,  1.5748031 , 0.5249344 , 1.0498688 , 2.3622047,  1.3123360,  1.3123360 , 1.0498688  ,1.5748031,  2.8871391 , 2.0997375,  1.8372703 , 3.1496063,  3.4120735 , 5.2493438 , 6.0367454, 11.0236220, 15.2230971, 18.6351706, 18.8976378),seq)

classIIIx <- data.frame(c( 0.7575758  ,0.7575758,  3.0303030,  4.5454545 , 5.3030303 , 5.3030303, 14.3939394 , 5.3030303 ,11.3636364 , 5.3030303 ,12.1212121, 11.3636364 , 3.7878788 , 0.7575758 , 1.5151515,  1.5151515, 0, 3.0303030 , 0.7575758 , 9.0909091),seq)

classIx <- data.frame(c(5.6910569 , 8.1300813 , 5.6910569 , 1.6260163 , 8.1300813 , 9.7560976 , 4.8780488 , 7.3170732 , 7.3170732  ,4.8780488 , 8.9430894,  4.8780488,0,  3.2520325,0 , 0.8130081 , 1.6260163,  3.2520325 ,3.2520325 ,10.5691057),seq)

classIIx <- data.frame(c( 3.260870 , 4.347826,  4.347826  ,1.086957,  3.260870 , 5.434783 , 5.434783  ,3.260870,  7.608696 , 5.434783,  6.521739,0 , 2.173913,0,  4.347826,  3.260870 , 2.173913 , 2.173913,  4.347826  ,31.521739),seq)

ALG34x <- data.frame(c(2.255639, 3.007519, 3.007519, 2.255639, 4.511278, 3.759398, 9.022556 ,5.263158 ,6.015038, 9.022556, 5.263158, 4.511278, 5.2631588, 9.022556, 5.263158 ,6.015038, 6.015038 ,4.511278, 0,0),seq)



all <- data.frame(classIx$seq, classIx$c.5.6910569..8.1300813..5.6910569..1.6260163..8.1300813..9.7560976.., classIIx$c.3.26087..4.347826..4.347826..1.086957..3.26087..5.434783..5.434783..  , classIIIx$c.0.7575758..0.7575758..3.030303..4.5454545..5.3030303..5.3030303..)
colnames(all) <- c("label", "Class-I", "Class-II", "Class-III")
print(all)

all2 <- melt(all)

p <- ggplot(all2, aes(x=label, y=value, fill=variable)) + geom_bar(position="dodge", stat = "identity") + scale_fill_manual(values=c("Grey", "Black","Blue")) + theme_classic() + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins log10 (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets_ChIPseqClasses.pdf", p, height = 10, width = 15, useDingbats=FALSE)



piRNAtop1 <- subset(piRNAtop, piRNAtop$V1 %in% classI$V1)
piRNAtop2 <- subset(piRNAtop, piRNAtop$V1 %in% classII$V1)
piRNAtop3 <- subset(piRNAtop, piRNAtop$V1 %in% classIII$V1)
piRNAtop4 <- subset(piRNAtop, piRNAtop$V1 %in% up$genes)


ALG34top1 <- subset(ALG34, ALG34$V1 %in% classI$V1)
ALG34top2 <- subset(ALG34, ALG34$V1 %in% classII$V1)
ALG34top3 <- subset(ALG34, ALG34$V1 %in% classIII$V1)


CSR1top1 <-  subset(CSR1, CSR1$V1 %in% classI$V1)
CSR1top2 <-  subset(CSR1, CSR1$V1 %in% classII$V1)
CSR1top3 <-  subset(CSR1, CSR1$V1 %in% classIII$V1)





#load the normalized data

data <- read.csv("RPB9ribo_mRNA_smallRNA_bins_log10TPM.csv", header = TRUE, stringsAsFactors = FALSE)

head(data)


# first, subset with piRNA targets and other arganoute targets


piRNAtop <- read.delim("./gene_targets/piRNAtargets.txt", header = FALSE)
piRNA2MM<- read.delim("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs/2MM_piRNatargets.txt", header = TRUE)
CSR1 <- read.delim("./gene_targets/CSR-1.txt", header = FALSE)
HRDE1 <- read.delim("./gene_targets/HRDE1.txt", header = FALSE)
ERGO1 <- read.delim("./gene_targets/ERGO-1.txt", header = FALSE)
ALG34 <- read.delim("./gene_targets/ALG-34.txt", header = FALSE)
germlineRNAi <- read.delim("./gene_targets/germline.txt", header = FALSE)
somaRNAi <- read.delim("./gene_targets/soma.txt", header = FALSE)
wago <- read.delim("./gene_targets/WAGO.txt", header = FALSE)


data_wago <- subset(data, data$genes %in% wago$V1)
freqwago <- as.data.frame(table(data_wago$label))
freqwago$abundance <- (freqwago$Freq/sum(freqwago$Freq))*100
write.csv(data_wago, "ribodata_wago.csv")

data_piRNA <- subset(data, data$genes %in% piRNAtop$V1)
freqpiRNA <- as.data.frame(table(data_piRNA$label))
freqpiRNA$abundance <- (freqpiRNA$Freq/sum(freqpiRNA$Freq))*100
write.csv(freqpiRNA, "./gene_targets/Freq_piRNA_Ribo.csv")
write.csv(data_piRNA, "ribodata_piRNA.csv")

data_piRNA2MM <- subset(data, data$genes %in% piRNA2MM$Gene.stable.ID)
freqpiRNA2MM <- as.data.frame(table(data_piRNA2MM$label))
freqpiRNA2MM$abundance <- (freqpiRNA2MM$Freq/sum(freqpiRNA2MM$Freq))*100
write.csv(data_piRNA2MM, "ribodata_piRNA2MM.csv")

data_CSR1 <- subset(data, data$genes %in% CSR1$V1)
freqCSR1 <- as.data.frame(table(data_CSR1$label))
freqCSR1$abundance <- (freqCSR1$Freq/sum(freqCSR1$Freq))*100
write.csv(freqCSR1, "./gene_targets/Freq_CSR1_Ribo.csv")
write.csv(data_CSR1, "ribodata_CSR1.csv")

data_HRDE1 <- subset(data, data$genes %in% HRDE1$V1)
freqHRDE1 <- as.data.frame(table(data_HRDE1$label))
freqHRDE1$abundance <- (freqHRDE1$Freq/sum(freqHRDE1$Freq))*100
write.csv(freqHRDE1, "./gene_targets/Freq_HRDE1_Ribo.csv")
write.csv(data_HRDE1, "ribodata_HRDE1.csv")

data_ERGO1 <- subset(data, data$genes %in% ERGO1$V1)
freqERGO1 <- as.data.frame(table(data_ERGO1$label))
freqERGO1$abundance <- (freqERGO1$Freq/sum(freqERGO1$Freq))*100
write.csv(data_ERGO1, "ribodata_ERGO1.csv")


data_ALG34 <- subset(data, data$genes %in% ALG34$V1)
freqALG34 <- as.data.frame(table(data_ALG34$label))
freqALG34$abundance <- (freqALG34$Freq/sum(freqALG34$Freq))*100
write.csv(data_ALG34, "ribodata_ALG34.csv")


data_germlineRNAi <- subset(data, data$genes %in% germlineRNAi$V1)
freqgermlineRNAi <- as.data.frame(table(data_germlineRNAi$label))
freqgermlineRNAi$abundance <- (freqgermlineRNAi$Freq/sum(freqgermlineRNAi$Freq))*100
write.csv(data_germlineRNAi, "ribodata_germlineRNAi.csv")


data_somaRNAi <- subset(data, data$genes %in% somaRNAi$V1)
freqsomaRNAi <- as.data.frame(table(data_somaRNAi$label))
freqsomaRNAi$abundance <- (freqsomaRNAi$Freq/sum(freqsomaRNAi$Freq))*100
write.csv(data_somaRNAi, "ribodata_somaRNAi.csv")


seq <- as.vector(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U"))

piRNAx <- data.frame(c(0.3546099,0 , 0.3546099 , 0.3546099 , 1.0638298 , 1.4184397 , 0.7092199 , 3.5460993 , 1.4184397,  1.4184397,  3.1914894 , 1.7730496 , 1.7730496 , 1.4184397 , 4.2553191 , 5.6737589  ,4.9645390 ,10.6382979 ,18.7943262, 36.8794326),seq)

CSR1x <- data.frame(c( 0.2597403,0,0,  0.2597403,0,  0.2597403 , 0.5194805 , 1.5584416 , 4.9350649 , 5.9740260 , 7.2727273 , 7.2727273 , 9.3506494 , 9.8701299 , 9.6103896 ,10.1298701  ,9.8701299 , 6.7532468 ,10.9090909 ,5.1948052),seq)

HRDE1x <- data.frame(c(0,0,0,0,0,0,0.2793296,0,0,  0.2793296 , 0.5586592,  0.2793296,0,  0.5586592,  0.2793296,  1.1173184,  1.6759777,  4.1899441, 15.3631285, 75.4189944),seq)

ERGO1x <- data.frame(c(0,0,0,0,0,0,0,0,0,0,0, 0, 0  ,0,0,0 , 4.545455,  4.545455,  4.545455, 86.363636),seq)

wagox <- data.frame(c(0.1199041,0,  0.1199041 ,0, 0.2398082 , 0.3597122,  0.5995204,  1.1990408 , 1.4388489 , 1.5587530  ,3.2374101  ,3.3573141 , 3.5971223  ,3.1175060 , 3.9568345 , 5.1558753 , 6.9544365,  9.4724221, 16.1870504 ,39.3285372),seq)

ALG34x <- data.frame(c(0,0,0,1.357466, 2.262443, 4.072398, 5.882353, 4.977376, 6.787330, 9.502262, 6.334842, 5.882353, 6.334842, 5.882353 ,8.144796, 7.692308, 3.619910 ,9.049774, 6.787330, 5.429864),seq)


all <- data.frame(piRNAx$seq, piRNAx$c.0.3546099..0..0.3546099..0.3546099..1.0638298..1.4184397..0.7092199.., wagox$c.0.1199041..0..0.1199041..0..0.2398082..0.3597122..0.5995204.., HRDE1x$c.0..0..0..0..0..0..0.2793296..0..0..0.2793296..0.5586592..0.2793296.., ERGO1x$c.0..0..0..0..0..0..0..0..0..0..0..0..0..0..0..0..4.545455..4.545455.., ALG34x$c.0..0..0..1.357466..2.262443..4.072398..5.882353..4.977376.., CSR1x$c.0.2597403..0..0..0.2597403..0..0.2597403..0.5194805..1.5584416..)
colnames(all) <- c("label", "piRNA", "WAGO", "HRDE-1", "ERGO1","ALG3/4","CSR1")
print(all)

all2 <- melt(all)

p <- ggplot(all2, aes(x=label, y=value, fill=variable)) + geom_bar(position="dodge", stat = "identity") + theme_classic() + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins log10 (increasing -->)")

ggsave("./gene_targets/ribo_mj261_mRNAsmallRNA_genetargets.pdf", p, height = 10, width = 15, useDingbats=FALSE)


all <- data.frame(piRNAx$seq, piRNAx$c.0.3546099..0..0.3546099..0.3546099..1.0638298..1.4184397..0.7092199.., HRDE1x$c.0..0..0..0..0..0..0.2793296..0..0..0.2793296..0.5586592..0.2793296..,  CSR1x$c.0.2597403..0..0..0.2597403..0..0.2597403..0.5194805..1.5584416..)
colnames(all) <- c("label", "piRNA", "HRDE-1","CSR1")
print(all)

all2 <- melt(all)

p <- ggplot(all2, aes(x=label, y=value, fill=variable)) + geom_bar(position="dodge", stat = "identity") + theme_classic() + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins log10 (increasing -->)")

ggsave("./gene_targets/ribo_mj261_mRNAsmallRNA_genetargets_onlyHRDE1CSR1.pdf", p, height = 10, width = 15, useDingbats=FALSE)


#load the normalized data

data <- read.csv("RPB9ribo_mRNA_smallRNA_bins_log10TPM.csv", header = TRUE, stringsAsFactors = FALSE)

head(data)


# then, re-cluster the last 4 bins. this might give a more sensitive information. 


piRNAtop <- read.delim("./gene_targets/piRNAtargets.txt", header = FALSE)
CSR1 <- read.delim("./gene_targets/CSR-1.txt", header = FALSE)
ERGO1 <- read.delim("./gene_targets/ERGO-1.txt", header = FALSE)
ALG34 <- read.delim("./gene_targets/ALG-34.txt", header = FALSE)
germlineRNAi <- read.delim("./gene_targets/germline.txt", header = FALSE)
somaRNAi <- read.delim("./gene_targets/soma.txt", header = FALSE)
wago <- read.delim("./gene_targets/WAGO.txt", header = FALSE)


data_wago <- subset(data, data$genes %in% wago$V1)
freqwago <- as.data.frame(table(data_wago$label))
freqwago$abundance <- (freqwago$Freq/sum(freqwago$Freq))*100
write.csv(data_wago, "data_wagosub.csv")

data_piRNA <- subset(data, data$genes %in% piRNAtop$V1)
freqpiRNA <- as.data.frame(table(data_piRNA$label))
freqpiRNA$abundance <- (freqpiRNA$Freq/sum(freqpiRNA$Freq))*100
write.csv(data_piRNA, "data_piRNAsub.csv")

data_CSR1 <- subset(data, data$genes %in% CSR1$V1)
freqCSR1 <- as.data.frame(table(data_CSR1$label))
freqCSR1$abundance <- (freqCSR1$Freq/sum(freqCSR1$Freq))*100
write.csv(data_CSR1, "data_CSR1sub.csv")


data_ERGO1 <- subset(data, data$genes %in% ERGO1$V1)
freqERGO1 <- as.data.frame(table(data_ERGO1$label))
freqERGO1$abundance <- (freqERGO1$Freq/sum(freqERGO1$Freq))*100
write.csv(data_ERGO1, "data_ERGO1sub.csv")


data_ALG34 <- subset(data, data$genes %in% ALG34$V1)
freqALG34 <- as.data.frame(table(data_ALG34$label))
freqALG34$abundance <- (freqALG34$Freq/sum(freqALG34$Freq))*100
write.csv(data_ALG34, "data_ALG34sub.csv")


data_germlineRNAi <- subset(data, data$genes %in% germlineRNAi$V1)
freqgermlineRNAi <- as.data.frame(table(data_germlineRNAi$label))
freqgermlineRNAi$abundance <- (freqgermlineRNAi$Freq/sum(freqgermlineRNAi$Freq))*100
write.csv(data_germlineRNAi, "data_germlineRNAisub.csv")


data_somaRNAi <- subset(data, data$genes %in% somaRNAi$V1)
freqsomaRNAi <- as.data.frame(table(data_somaRNAi$label))
freqsomaRNAi$abundance <- (freqsomaRNAi$Freq/sum(freqsomaRNAi$Freq))*100
write.csv(data_somaRNAi, "data_somaRNAisub.csv")


seq <- as.vector(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U"))

piRNAx <- data.frame(c(0.7874016,  1.5748031 , 0.5249344 , 1.0498688 , 2.3622047,  1.3123360,  1.3123360 , 1.0498688  ,1.5748031,  2.8871391 , 2.0997375,  1.8372703 , 3.1496063,  3.4120735 , 5.2493438 , 6.0367454, 11.0236220, 15.2230971, 18.6351706, 18.8976378),seq)

CSR1x <- data.frame(c(0,0.2631579 , 0.7894737 , 2.3684211 , 1.8421053 , 3.4210526,  3.4210526 , 5.0000000 , 7.8947368,  7.8947368  ,7.6315789 , 6.3157895, 10.5263158,  8.6842105 ,11.0526316 , 7.6315789,  7.8947368  ,5.5263158 , 1.8421053,0),seq)

ERGO1x <- data.frame(c(0,0,1.960784,0,3.921569,0,0,0,0,0,1.960784,3.921569,0,1.960784,  1.960784 , 3.921569,  3.921569 ,17.647059, 19.607843, 39.215686),seq)

wagox <- data.frame(c(0,0,0,0, 0.09624639 ,0, 0.09624639 ,0, 0.09624639 ,0, 0.38498556,  0.48123195,  1.25120308  ,2.21366699  ,2.98363811,  5.00481232  ,5.38979788,  9.43214629 ,20.01924928 ,52.55052936),seq)

ALG34x <- data.frame(c(0.2518892 , 0.2518892,  0.2518892 , 0.2518892 , 3.5264484 , 4.0302267 , 6.2972292,  9.5717884 , 8.3123426 , 8.5642317 ,11.8387909 ,10.3274559,  5.5415617 , 2.2670025  ,6.5491184  ,4.0302267, 4.0302267 , 5.2896725  ,4.7858942 , 4.0302267),seq)


all <- data.frame(piRNAx$seq, piRNAx$c.0..0..0..0..0.2645503..0..0..0.2645503..0.7936508..0.5291005.., wagox$c.0..0..0..0..0.09624639..0..0.09624639..0..0.09624639..0..0.38498556..  , HRDE1x$c.0..0..0..0..0..0..0..0..0.2197802..0..0..0.2197802..0.6593407.., ERGO1x$c.0..0..0..0..0..0..0..0..0..0..0..2.439024..2.439024..2.439024.., ALG34x$c.0.2518892..0.2518892..0.2518892..0.2518892..3.5264484..4.0302267..,CSR1x$c.0.2590674..0..0..0..0..0.5181347..0..0..0..0..0..0.7772021..)
colnames(all) <- c("label", "piRNA", "WAGO", "HRDE-1", "ERGO-1","ALG-3/4","CSR-1")
print(all)

all2 <- melt(all)

p <- ggplot(all2, aes(x=label, y=value, fill=variable)) + geom_bar(position="dodge", stat = "identity") + theme_classic() + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins log10 (increasing -->)")

ggsave("./gene_targets/polyA_mj261_mRNAsmallRNA_genetargets.pdf", p, height = 10, width = 15, useDingbats=FALSE)

# then, repeat this for upregulated and ChIP-seq CLassI, II, III genes

polyA <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_ribozero/mj261_ribozero_diff_osc.csv", header = TRUE, stringsAsFactors = FALSE)

up <- polyA[polyA$threshold =="A", ]

up <- na.omit(up)

down <- polyA[polyA$threshold =="B", ]

down <- na.omit(down)


data <- read.csv("RPB9ribo_mRNA_smallRNA_bins_log10TPM.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)

data_cluster <- read.csv("RPB9ribo_mRNA_smallRNA_bins_log10TPM_last4bins.csv", header = TRUE, stringsAsFactors = FALSE)
head(data_cluster)


classI <- read.delim("upregulated_classI_ribo.tsv", header = FALSE)
classII <- read.delim("upregulated_classII_ribo.tsv", header = FALSE)
classIII <- read.delim("upregulated_classIII_ribo.tsv", header = FALSE)

data_up <- subset(data, data$genes %in% up$genes)
data_down <- subset(data, data$genes %in% down$genes)
data_upcluster <- subset(data_cluster, data_cluster$genes %in% up$genes)
data_downcluster <- subset(data_cluster, data_cluster$genes %in% down$genes)
data_classI <- subset(data, data$genes %in% classI$V1)
data_classII <- subset(data, data$genes %in% classII$V1)
data_classIII <- subset(data, data$genes %in% classIII$V1)

data_classI2 <- subset(data_cluster, data_cluster$genes %in% classI$V1)
data_classII2 <- subset(data_cluster, data_cluster$genes %in% classII$V1)
data_classIII2 <- subset(data_cluster, data_cluster$genes %in% classIII$V1)


freqdown <- as.data.frame(table(data_down$label))
freqdown$abundance <- (freqdown$Freq/sum(freqdown$Freq))*100
write.csv(data_down, "ribodata_down.csv")

p <- ggplot(freqdown, aes(x=Var1, y=abundance)) + geom_bar( stat = "identity", color = "Black",width =0.25) + theme_classic() + scale_x_discrete(limits = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U")) + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins (increasing -->)")

ggsave("./gene_targets/ribo_mj261_mRNAsmallRNA_genetargets_downregulatedgenes.pdf", p, height = 10, width = 15, useDingbats=FALSE)


freqdown2 <- as.data.frame(table(data_downcluster$label))
freqdown2$abundance <- (freqdown2$Freq/sum(freqdown2$Freq))*100
write.csv(data_downcluster, "ribodata_down_cluster.csv")


p <- ggplot(freqdown2, aes(x=Var1, y=abundance)) + geom_bar( stat = "identity", color = "Black",width =0.25) + theme_classic() + scale_x_discrete(limits = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U")) + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins (increasing -->)")

ggsave("./gene_targets/ribo_mj261_mRNAsmallRNA_genetargets_downregulatedgenes_last4clusters.pdf", p, height = 10, width = 15, useDingbats=FALSE)



frequp <- as.data.frame(table(data_up$label))
frequp$abundance <- (frequp$Freq/sum(frequp$Freq))*100
write.csv(frequp, "./gene_targets/Freq_upregulated_Ribo.csv")
write.csv(data_up, "ribodata_up.csv")

p <- ggplot(frequp, aes(x=Var1, y=abundance)) + geom_bar( stat = "identity", color = "Black",width =0.25) + theme_classic() + scale_x_discrete(limits = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U")) + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins (increasing -->)")

ggsave("./gene_targets/ribo_mj261_mRNAsmallRNA_genetargets_upregulatedgenes.pdf", p, height = 10, width = 15, useDingbats=FALSE)


frequp2 <- as.data.frame(table(data_upcluster$label))
frequp2$abundance <- (frequp2$Freq/sum(frequp2$Freq))*100
write.csv(data_upcluster, "ribodata_up_cluster.csv")


p <- ggplot(frequp2, aes(x=Var1, y=abundance)) + geom_bar( stat = "identity", color = "Black",width =0.25) + theme_classic() + scale_x_discrete(limits = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U")) + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins (increasing -->)")

ggsave("./gene_targets/ribo_mj261_mRNAsmallRNA_genetargets_upregulatedgenes_last4clusters.pdf", p, height = 10, width = 15, useDingbats=FALSE)


freqclassI <- as.data.frame(table(data_classI$label))
freqclassI$abundance <- (freqclassI$Freq/sum(freqclassI$Freq))*100
write.csv(freqclassI, "./gene_targets/Freq_ClassIII_Ribo.csv")
write.csv(data_classI, "ribodata_classI.csv")

freqclassII <- as.data.frame(table(data_classII$label))
freqclassII$abundance <- (freqclassII$Freq/sum(freqclassII$Freq))*100
write.csv(freqclassII, "./gene_targets/Freq_ClassI_Ribo.csv")
write.csv(data_classII, "ribodata_classII.csv")


freqclassIII <- as.data.frame(table(data_classIII$label))
freqclassIII$abundance <- (freqclassIII$Freq/sum(freqclassIII$Freq))*100
write.csv(freqclassIII, "./gene_targets/Freq_ClassII_Ribo.csv")
write.csv(data_classIII, "ribodata_classIII.csv")

freqclassI2 <- as.data.frame(table(data_classI2$label))
freqclassI2$abundance <- (freqclassI2$Freq/sum(freqclassI2$Freq))*100
write.csv(data_classI2, "data_classIcluster.csv")

freqclassII2 <- as.data.frame(table(data_classII2$label))
freqclassII2$abundance <- (freqclassII2$Freq/sum(freqclassII2$Freq))*100
write.csv(data_classII2, "data_classIIcluster.csv")


freqclassIII2 <- as.data.frame(table(data_classIII2$label))
freqclassIII2$abundance <- (freqclassIII2$Freq/sum(freqclassIII2$Freq))*100
write.csv(data_classIII2, "data_classIIIcluster.csv")


seq <- as.vector(c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "R", "S", "T", "U"))

upx <- data.frame(c(0.7874016,  1.5748031 , 0.5249344 , 1.0498688 , 2.3622047,  1.3123360,  1.3123360 , 1.0498688  ,1.5748031,  2.8871391 , 2.0997375,  1.8372703 , 3.1496063,  3.4120735 , 5.2493438 , 6.0367454, 11.0236220, 15.2230971, 18.6351706, 18.8976378),seq)

classIx <- data.frame(c( 0, 5.952381 , 7.142857 , 4.761905,  7.142857  ,4.761905 , 7.142857, 11.904762 , 9.523810 , 8.333333, 10.714286 , 2.380952 , 3.571429,  3.571429 , 3.571429 , 2.380952 , 1.190476, 0 ,1.190476,4.761905),seq)

classIIx <- data.frame(c(5.494505 , 8.791209, 10.989011, 15.384615 , 7.692308,  9.890110,  7.692308,  4.395604,  4.395604 , 4.395604,0,  1.098901, 0, 1.098901 , 1.098901, 0, 1.098901,  2.197802  ,4.395604 , 9.890110),seq)

classIIIx <- data.frame(c(0,0,2.325581,0 , 4.651163 , 2.325581 , 4.651163,  6.976744,  2.325581,  4.651163,  2.325581,  2.325581,  2.325581, 0, 6.976744,  2.325581,  2.325581 ,11.627907, 11.627907, 30.232558),seq)

ALG34x <- data.frame(c(2.255639, 3.007519, 3.007519, 2.255639, 4.511278, 3.759398, 9.022556 ,5.263158 ,6.015038, 9.022556, 5.263158, 4.511278, 5.2631588, 9.022556, 5.263158 ,6.015038, 6.015038 ,4.511278, 0,0),seq)



all <- data.frame(classIx$seq, classIx$c.0..5.952381..7.142857..4.761905..7.142857..4.761905..7.142857.., classIIx$c.5.494505..8.791209..10.989011..15.384615..7.692308..9.89011..  , classIIIx$c.0..0..2.325581..0..4.651163..2.325581..4.651163..6.976744..)
colnames(all) <- c("label", "Class-I", "Class-II", "Class-III")
print(all)

all2 <- melt(all)

p <- ggplot(all2, aes(x=label, y=value, fill=variable)) + geom_bar(position="dodge", stat = "identity") + scale_fill_manual(values=c("Grey", "Black","Blue")) + theme_classic() + scale_y_continuous(limits = c(0,100)) + ylab("Percentage(%)") + xlab("small RNA density bins log10 (increasing -->)")

ggsave("./gene_targets/ribo_mj261_mRNAsmallRNA_genetargets_ChIPseqClasses.pdf", p, height = 10, width = 15, useDingbats=FALSE)


piRNAtop1 <- subset(piRNAtop, piRNAtop$V1 %in% classI$V1)
piRNAtop2 <- subset(piRNAtop, piRNAtop$V1 %in% classII$V1)
piRNAtop3 <- subset(piRNAtop, piRNAtop$V1 %in% classIII$V1)
piRNAtop4 <- subset(piRNAtop, piRNAtop$V1 %in% up$genes)


ALG34top1 <- subset(ALG34, ALG34$V1 %in% classI$V1)
ALG34top2 <- subset(ALG34, ALG34$V1 %in% classII$V1)
ALG34top3 <- subset(ALG34, ALG34$V1 %in% classIII$V1)


CSR1top1 <-  subset(CSR1, CSR1$V1 %in% classI$V1)
CSR1top2 <-  subset(CSR1, CSR1$V1 %in% classII$V1)
CSR1top3 <-  subset(CSR1, CSR1$V1 %in% classIII$V1)

