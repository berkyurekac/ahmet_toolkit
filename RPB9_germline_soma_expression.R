
library(dplyr)
library(reshape2)
library(ggplot2)
##### 
# RPB-9 germline vs. soma expression - from adult animals reinke 2018?
#####
setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_1")


# introduce total counts

counts <- read.csv("counts_germline_soma.csv", header = TRUE)
head(counts)


germ1 <- read.delim("GSM3268042_IGN_RNA-seq_rep1_expnew.txt", header = FALSE)
germ2 <- read.delim("GSM3268043_IGN_RNA-seq_rep2_expnew.txt", header = FALSE)
soma1 <- read.delim("GSM3268044_SOM_RNA-seq_rep1_expnew.txt", header = FALSE)
soma2 <- read.delim("GSM3268045_SOM_RNA-seq_rep2_expnew.txt", header = FALSE)

head(germ1)

frame <- data.frame(germ1[,c(1,2,6)], germ2[,6], soma1[,6], soma2[,6])
print(frame)
colnames(frame) <- c("genes", "length", "germ1", "germ2", "soma1", "soma2")
head(frame)

frame$germ1 <- (frame$germ1/sum(counts$X42_IGN_RNA.seq_rep1.bw))*1e9
frame$germ2 <- (frame$germ2/sum(counts$X43_IGN_RNA.seq_rep2.bw))*1e9
frame$soma1 <- (frame$soma1/sum(counts$X44_SOM_RNA.seq_rep1.bw))*1e9
frame$soma2 <- (frame$soma2/sum(counts$X45_SOM_RNA.seq_rep2.bw))*1e9

head(frame)
                        
a <- melt(frame)
a$value <- log2(a$value)

germline <- as.vector(c(5.889825,8.637745))
soma <- as.vector(c(6.639975, 4.421156))
x <- data.frame(germline, soma)

y <- melt(x)

p <- ggplot(y, aes(x=y$variable, y=y$value)) + geom_boxplot( color="Black", size=1) + theme_classic()  + scale_y_continuous(limits = c(0,10))

ggsave("RPB9_germline_vs_somaexpression_boxplot.pdf", p, height = 10, width = 10, useDingbats=FALSE)


germline <- as.vector(c(8.520147,9.5769540))
soma <- as.vector(c(0.224648, 0.992146))
x <- data.frame(germline, soma)

y <- melt(x)

p <- ggplot(y, aes(x=y$variable, y=y$value)) + geom_boxplot( color="Black", size=1) + theme_classic()  + scale_y_continuous(limits = c(0,10))

ggsave("HRDE1_germline_vs_somaexpression_boxplot.pdf", p, height = 10, width = 10, useDingbats=FALSE)


germline <- as.vector(c(7.40751744,8.54761807))
soma <- as.vector(c(8.04761807, 8.88259438))
x <- data.frame(germline, soma)

y <- melt(x)

p <- ggplot(y, aes(x=y$variable, y=y$value)) + geom_boxplot( color="Black", size=1) + theme_classic()  + scale_y_continuous(limits = c(0,10))

ggsave("TBA2_germline_vs_somaexpression_boxplot.pdf", p, height = 10, width = 10, useDingbats=FALSE)


germline <- as.vector(c(0.65111211,0.02491853))
soma <- as.vector(c(3.03676697, 3.86505988))
x <- data.frame(germline, soma)

y <- melt(x)

p <- ggplot(y, aes(x=y$variable, y=y$value)) + geom_boxplot( color="Black", size=1) + theme_classic()  + scale_y_continuous(limits = c(0,10))

ggsave("NRDE3_germline_vs_somaexpression_boxplot.pdf", p, height = 10, width = 10, useDingbats=FALSE)



##### 
# tissue specific expression from Ahringer and Serizay
#####
list <- read.delim("gene-list_full-report.txt", header = T)

head(list)

listrpb9 <- subset(list, list$Associated_gene_Name == "rpb-9")

listtba2 <- subset(list, list$Associated_gene_Name == "tba-2")


listgapdh <- subset(list, list$Associated_gene_Name == "gpd-4")

listhrde1 <- subset(list, list$Associated_gene_Name == "hrde-1")

listcdc42 <- subset(list, list$Associated_gene_Name == "cdc-42")

print(listrpb9)
print(listtba2)
print(listgapdh)
print(listhrde1)

rpb9 <- c(133.706 , 36.682 , 36.499 , 29.767 , 27.424)

tba2 <- c(693.294 , 686.725 , 417.916 , 343.08 , 338.958)

gapdh <- c(84.231,33.408,22.286,20.97,16.117)

hrde1 <- c(543.857 , 221.077 , 212.075 , 141.69 , 119.872)

cdc42 <- c(277.917 , 180.442 , 162.59 , 162.013 , 117.969)

names <- c("germline", "Intestine", "Hypodermis", "Muscle", "Neurons")


frame <- data.frame(rpb9, tba2, gapdh,hrde1, cdc42, names)

frame$norm <- frame$rpb9/frame$tba2

frame$norm2 <- frame$rpb9/frame$gapdh

frame$norm3 <- frame$rpb9/frame$cdc42

print(frame)




g <- ggplot(frame, aes(x=names, y=norm3)) + geom_bar(stat = "identity") + theme_classic() + scale_y_continuous(limits = c(0,0.6)) + ylab('rpb-9 TPM / cdc-42 TPM')

ggsave("rpb9_tissue_normalized_gapdh_new2.pdf", g, height = 5, width = 5, useDingbats=FALSE)


g <- ggplot(frame, aes(x=names, y=log10(rpb9))) + geom_bar(stat = "identity") + theme_classic()

ggsave("rpb9_tissue_rpb9.pdf", g, height = 5, width = 5, useDingbats=FALSE)


g <- ggplot(frame, aes(x=names, y=log10(gapdh))) + geom_bar(stat = "identity") + theme_classic()

ggsave("rpb9_tissue_gapdh.pdf", g, height = 5, width = 5, useDingbats=FALSE)


g <- ggplot(frame, aes(x=names, y=log10(tba2))) + geom_bar(stat = "identity") + theme_classic()

ggsave("rpb9_tissue_tba2.pdf", g, height = 5, width = 5, useDingbats=FALSE)










