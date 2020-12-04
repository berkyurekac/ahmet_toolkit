# this script was created by Ahmet Can Berkyurek on 06/04/2020, to see up and down genes on AMA-1 ChIP-seq!
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

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/others/subsets")

##### 
# introduce up and down for AMA-1 CHIP-seq
#####

# first for upregulated genes
upwt <- read.csv("AMA1_upregulated_wt.csv", header = TRUE, stringsAsFactors = FALSE)
head(upwt)
upmj261 <- read.csv("AMA1_upregulated_mj261.csv", header = TRUE, stringsAsFactors = FALSE)
head(upmj261)

upmerged <- merge(upwt, upmj261, by.x="metadata_name", by.y="metadata_name")
head(upmerged)

upclassI <- subset(upmerged, upmerged$RowMeans.x > 0  & upmerged$RowMeans.x > upmerged$RowMeans.y) # RNA pol II higher in wt than in mj261. all values > 0
write.csv(upclassI, "AMA1ChIPseq_upregulatedgenes_wt_bigger_mj261_CLASS1.csv")

upclassII <- subset(upmerged, upmerged$RowMeans.y > 0  & upmerged$RowMeans.y > (upmerged$RowMeans.x)*1.5) # RNA pol II higher in mj261 than in wt. all values > 0
write.csv(upclassII, "AMA1ChIPseq_upregulatedgenes_mj261_bigger_wt_CLASS2.csv")

a <- subset(upmerged, !(upmerged$metadata_name) %in% upclassII$metadata_name)

b <- subset(a, !(a$metadata_name) %in% upclassI$metadata_name) 

upclassIII <- subset(b, !(b$RowMeans.y) > (b$RowMeans.x)*1.1) # RNA pol II relatively similar in wt and mj261. all values > 0
write.csv(upclassIII, "AMA1ChIPseq_upregulatedgenes_mj261_equal_wt_CLASS3.csv")


# first for downregulated genes
downwt <- read.csv("AMA1_downregulated_wt.csv", header = TRUE, stringsAsFactors = FALSE)
head(downwt)
downmj261 <- read.csv("AMA1_downregulated_mj261.csv", header = TRUE, stringsAsFactors = FALSE)
head(downmj261)

downmerged <- merge(downwt, downmj261, by.x="metadata_name", by.y="metadata_name")
head(downmerged)

downclassI <- subset(downmerged, downmerged$RowMeans.x > 0  & downmerged$RowMeans.x > downmerged$RowMeans.y) # RNA pol II higher in wt than in mj261. all values > 0
write.csv(downclassI, "AMA1ChIPseq_downregulatedgenes_wt_bigger_mj261_CLASS1.csv")

downclassII <- subset(downmerged, downmerged$RowMeans.y > 0  & downmerged$RowMeans.y > (downmerged$RowMeans.x)*1.5) # RNA pol II higher in mj261 than in wt. all values > 0
write.csv(downclassII, "AMA1ChIPseq_downregulatedgenes_mj261_bigger_wt_CLASS2.csv")

a <- subset(downmerged, !(downmerged$metadata_name) %in% downclassII$metadata_name)

b <- subset(a, !(a$metadata_name) %in% downclassI$metadata_name) 

downclassIII <- subset(b, b$RowMeans.x >0 & b$RowMeans.y > 0) # RNA pol II relatively similar in wt and mj261. all values > 0
write.csv(downclassIII, "AMA1ChIPseq_downregulatedgenes_mj261_equal_wt_CLASS3.csv")


# subset diff and non-diff genes

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/All_genes")

allgenes <- read.delim("Genes_ps562_source_ce11_ranges_ANN001_Ae970951fixed3.bed", header = FALSE, stringsAsFactors = FALSE)

polyA <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA//mj261_polyA_diff_osc.csv", header = TRUE, stringsAsFactors = FALSE)
head(polyA)

diff <- subset(polyA, !(polyA$threshold) == "C")
head(diff)

diff_genome <- subset(allgenes, allgenes$V4 %in% diff$X)
head(diff_genome)

write.table(x = diff_genome, 
            file = "diff_genome.bed", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

non_diff <- subset(allgenes, !(allgenes$V4) %in% diff_genome$V4)
head(non_diff)

write.table(x = non_diff, 
            file = "non_diff_genome.bed", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')





setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/polyA_allreplicates")


##### 
# inroduce up and down for AMA-1 CHIP-seq
#####

# first for upregulated genes
upwt <- read.csv("AMA1_upregulated_wt_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)
head(upwt)
upmj261 <- read.csv("AMA1_upregulated_mj261_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)
head(upmj261)

upmerged <- merge(upwt, upmj261, by.x="metadata_name", by.y="metadata_name")
head(upmerged)

upclassI <- subset(upmerged, upmerged$RowMeans.x > 1  & upmerged$RowMeans.x > (upmerged$RowMeans.y)) # RNA pol II higher in wt than in mj261. all values > 0
write.csv(upclassI, "AMA1ChIPseq_upregulatedgenes_wt_bigger_mj261_CLASS1_allreplicates_new.csv")

upclassII <- subset(upmerged, upmerged$RowMeans.y > 0  & upmerged$RowMeans.y > (upmerged$RowMeans.x)*1.5) # RNA pol II higher in mj261 than in wt. all values > 0
write.csv(upclassII, "AMA1ChIPseq_upregulatedgenes_mj261_bigger_wt_CLASS2_allreplicates.csv")

a <- subset(upmerged, !(upmerged$metadata_name) %in% upclassII$metadata_name)

b <- subset(a, !(a$metadata_name) %in% upclassI$metadata_name) 

upclassIII <- subset(b, b$RowMeans.x >0 & b$RowMeans.y > 0) # RNA pol II relatively similar in wt and mj261. all values > 0
write.csv(upclassIII, "AMA1ChIPseq_upregulatedgenes_mj261_equal_wt_CLASS3_allreplicates.csv")


# first for downregulated genes
downwt <- read.csv("AMA1_downregulated_wt_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)
head(downwt)
downmj261 <- read.csv("AMA1_downregulated_mj261_allreplicates.csv", header = TRUE, stringsAsFactors = FALSE)
head(downmj261)

downmerged <- merge(downwt, downmj261, by.x="metadata_name", by.y="metadata_name")
head(downmerged)

downclassI <- subset(downmerged, downmerged$RowMeans.x > 0  & downmerged$RowMeans.x > downmerged$RowMeans.y) # RNA pol II higher in wt than in mj261. all values > 0
write.csv(downclassI, "AMA1ChIPseq_downregulatedgenes_wt_bigger_mj261_CLASS1_allreplicates.csv")

downclassII <- subset(downmerged, downmerged$RowMeans.y > 0  & downmerged$RowMeans.y > (downmerged$RowMeans.x)*1.5) # RNA pol II higher in mj261 than in wt. all values > 0
write.csv(downclassII, "AMA1ChIPseq_downregulatedgenes_mj261_bigger_wt_CLASS2_allreplicates.csv")

a <- subset(downmerged, !(downmerged$metadata_name) %in% downclassII$metadata_name)

b <- subset(a, !(a$metadata_name) %in% downclassI$metadata_name) 

downclassIII <- subset(b, b$RowMeans.x >0 & b$RowMeans.y > 0) # RNA pol II relatively similar in wt and mj261. all values > 0
write.csv(downclassIII, "AMA1ChIPseq_downregulatedgenes_mj261_equal_wt_CLASS3_allreplicates.csv")


# subset diff and non-diff genes

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/All_genes")

allgenes <- read.delim("Genes_ps562_source_ce11_ranges_ANN001_Ae970951fixed3.bed", header = FALSE, stringsAsFactors = FALSE)

polyA <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA//mj261_polyA_diff_osc.csv", header = TRUE, stringsAsFactors = FALSE)
head(polyA)

diff <- subset(polyA, !(polyA$threshold) == "C")
head(diff)

diff_genome <- subset(allgenes, allgenes$V4 %in% diff$X)
head(diff_genome)

write.table(x = diff_genome, 
            file = "diff_genome.bed", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

non_diff <- subset(allgenes, !(allgenes$V4) %in% diff_genome$V4)
head(non_diff)

write.table(x = non_diff, 
            file = "non_diff_genome.bed", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')


setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/subsets/ribozero")


##### 
# inroduce up and down for AMA-1 CHIP-seq
#####

# first for upregulated genes
upwt <- read.csv("AMA1_upregulated_wt_ribo.csv", header = TRUE, stringsAsFactors = FALSE)
head(upwt)
upmj261 <- read.csv("AMA1_upregulated_mj261_ribo.csv", header = TRUE, stringsAsFactors = FALSE)
head(upmj261)

upmerged <- merge(upwt, upmj261, by.x="metadata_name", by.y="metadata_name")
head(upmerged)

upclassI <- subset(upmerged, upmerged$RowMeans.x > 0  & upmerged$RowMeans.x > upmerged$RowMeans.y) # RNA pol II higher in wt than in mj261. all values > 0
write.csv(upclassI, "AMA1ChIPseq_upregulatedgenes_wt_bigger_mj261_CLASS1_ribo.csv")

upclassII <- subset(upmerged, upmerged$RowMeans.y > 0  & upmerged$RowMeans.y > upmerged$RowMeans.x) # RNA pol II higher in mj261 than in wt. all values > 0
write.csv(upclassII, "AMA1ChIPseq_upregulatedgenes_mj261_bigger_wt_CLASS2_ribo.csv")

a <- subset(upmerged, !(upmerged$metadata_name) %in% upclassII$metadata_name)

b <- subset(a, !(a$metadata_name) %in% upclassI$metadata_name) 

upclassIII <- subset(b, b$RowMeans.x >0 & b$RowMeans.y > 0) # RNA pol II relatively similar in wt and mj261. all values > 0
write.csv(upclassIII, "AMA1ChIPseq_upregulatedgenes_mj261_equal_wt_CLASS3_ribo.csv")


# first for downregulated genes
downwt <- read.csv("AMA1_downregulated_wt_ribo.csv", header = TRUE, stringsAsFactors = FALSE)
head(downwt)
downmj261 <- read.csv("AMA1_downregulated_mj261_ribo.csv", header = TRUE, stringsAsFactors = FALSE)
head(downmj261)

downmerged <- merge(downwt, downmj261, by.x="metadata_name", by.y="metadata_name")
head(downmerged)

downclassI <- subset(downmerged, downmerged$RowMeans.x > 0  & downmerged$RowMeans.x > downmerged$RowMeans.y) # RNA pol II higher in wt than in mj261. all values > 0
write.csv(downclassI, "AMA1ChIPseq_downregulatedgenes_wt_bigger_mj261_CLASS1_ribo.csv")

downclassII <- subset(downmerged, downmerged$RowMeans.y > 0  & downmerged$RowMeans.y > (downmerged$RowMeans.x)*1.5) # RNA pol II higher in mj261 than in wt. all values > 0
write.csv(downclassII, "AMA1ChIPseq_downregulatedgenes_mj261_bigger_wt_CLASS2_ribo.csv")

a <- subset(downmerged, !(downmerged$metadata_name) %in% downclassII$metadata_name)

b <- subset(a, !(a$metadata_name) %in% downclassI$metadata_name) 

downclassIII <- subset(b, b$RowMeans.x >0 & b$RowMeans.y > 0) # RNA pol II relatively similar in wt and mj261. all values > 0
write.csv(downclassIII, "AMA1ChIPseq_downregulatedgenes_mj261_equal_wt_CLASS3_ribo.csv")


# subset diff and non-diff genes

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/Figure_4/All_genes")

allgenes <- read.delim("Genes_ps562_source_ce11_ranges_ANN001_Ae970951fixed3.bed", header = FALSE, stringsAsFactors = FALSE)

polyA <- read.csv("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/RNA-seq_polyA//mj261_polyA_diff_osc.csv", header = TRUE, stringsAsFactors = FALSE)
head(polyA)

diff <- subset(polyA, !(polyA$threshold) == "C")
head(diff)

diff_genome <- subset(allgenes, allgenes$V4 %in% diff$X)
head(diff_genome)

write.table(x = diff_genome, 
            file = "diff_genome.bed", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

non_diff <- subset(allgenes, !(allgenes$V4) %in% diff_genome$V4)
head(non_diff)

write.table(x = non_diff, 
            file = "non_diff_genome.bed", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

