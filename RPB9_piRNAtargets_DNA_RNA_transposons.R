# this script was created by Ahmet Can Berkyurek on 16/09/2020 to compare the ratios of DNA and RNA transposons in piRNA targets!
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

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Final Figures/EMBO/revision/bioinformatics/piRNA_targets")


#####
# introduce all transposons and matching color codes, first
#####


TE <- read.delim("Repeats_ps562_source_ce11_ranges_ANN003_A9970983.bed", header = F, stringsAsFactors = F)

head(TE)

TE2 <- TE[!duplicated(TE[,c('V4','V9')]),]

freq <- as.data.frame(table(TE2$V9))
freq$abundance <- (freq$Freq/sum(freq$Freq))*100

print(freq)

# all transposons
other = 32
DNA = 98
RNA = 57



#####
# introduce piRNA target transposons, and find frequencies
#####

piRNAtargets <- read.delim("1220952s_Additional_data_file_5.txt", header = T, stringsAsFactors = F)
head(piRNAtargets)

piRNAtargets_type <- subset(TE2, TE2$V4 %in% piRNAtargets$Name)
head(piRNAtargets_type)


freq <- as.data.frame(table(piRNAtargets_type$V9))
freq$abundance <- (freq$Freq/sum(freq$Freq))*100

print(freq)

# piRNA target transposons
other = 10
DNA = 62
RNA = 30













