library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(AnnotationHub)
ahub <- AnnotationHub()
table(ahub$rdataclass)


##########
setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/ChIP-seq/Julie_Ahringer_data")

########## write.table(x = bdg, 
##########           file = "SX1316_rep1_log2.bdg", 
##########            row.names = F, 
##########            col.names = F, 
##########            quote = F, 
##########            sep = '\t')
##########
##########
##########

##log2 FE
chrsize <- read.delim("ce11.chrom.sizes", header = FALSE, stringsAsFactors = FALSE)
colnames(chrsize) <- c("seqnames", "seqlenghts")
print(chrsize)


file <- import.bedGraph("SRR4319294_chromatin_peaks_linear.bdg")
head(file)

file$score <- log2(file$score+0.1)
head(file)


fileseqinfo <- Seqinfo(chrsize$seqnames, chrsize$seqlenghts, isCircular=NA, genome="ce11")
seqlevels(file) <- seqlevels(fileseqinfo)
seqlevels(file)

seqinfo(file) <- fileseqinfo

export.bw(file, "SRR4319294_rep1_log2.bw")


setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/Hiseq/ollas-rpb9/new_normalisation/substraction")

N2_rep1 <- read.table("N2_rep1._coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(N2_rep1)

N2_rep2 <- read.table("N2_rep2_coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(N2_rep2)

N2_rep3 <- read.table("N2_rep3_coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(N2_rep3)

N2_rep1$mean <- (N2_rep1$V3+N2_rep2$V3+N2_rep3$V3)/3


SX3220_rep1 <- read.table("SX3220-rep1_coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(SX3220_rep1)
SX3220_rep1$V3 <- SX3220_rep1$V3-N2_rep1$mean
SX3220_rep1$end <- SX3220_rep1$V2+1
SX3220_rep1 <- SX3220_rep1[,c(1,2,4,3)]

write.table(x = SX3220_rep1, 
          file = "SX3220_rep1_substracted.bdg", 
           row.names = F, 
           col.names = F, 
          quote = F, 
           sep = '\t')

SX3220_rep2 <- read.table("SX3220-rep2_coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(SX3220_rep2)
SX3220_rep2$V3 <- SX3220_rep2$V3-N2_rep1$mean
SX3220_rep2$end <- SX3220_rep2$V2+1
SX3220_rep2 <- SX3220_rep2[,c(1,2,4,3)]

write.table(x = SX3220_rep2, 
            file = "SX3220_rep2_substracted.bdg", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

SX3220_rep3 <- read.table("SX3220_rep3_coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(SX3220_rep3)
SX3220_rep3$V3 <- SX3220_rep3$V3-N2_rep1$mean
SX3220_rep3$end <- SX3220_rep3$V2+1
SX3220_rep3 <- SX3220_rep3[,c(1,2,4,3)]

write.table(x = SX3220_rep3, 
            file = "SX3220_rep3_substracted.bdg", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

SX3221_rep1 <- read.table("SX3221_rep1_coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(SX3221_rep1)
SX3221_rep1$V3 <- SX3221_rep1$V3-N2_rep1$mean
SX3221_rep1$end <- SX3221_rep1$V2+1
SX3221_rep1 <- SX3221_rep1[,c(1,2,4,3)]

write.table(x = SX3221_rep1, 
            file = "SX3221_rep1_substracted.bdg", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

SX3221_rep2 <- read.table("SX3221_rep2_coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(SX3221_rep2)
SX3221_rep2$V3 <- SX3221_rep2$V3-N2_rep1$mean
SX3221_rep2$end <- SX3221_rep2$V2+1
SX3221_rep2 <- SX3221_rep2[,c(1,2,4,3)]

write.table(x = SX3221_rep2, 
            file = "SX3221_rep2_substracted.bdg", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

SX3221_rep3 <- read.table("SX3221_rep3_coverage_raw_2.bdg", header = FALSE, stringsAsFactors = FALSE)
head(SX3221_rep3)
SX3221_rep3$V3 <- SX3221_rep3$V3-N2_rep1$mean
SX3221_rep3$end <- SX3221_rep3$V2+1
SX3221_rep3 <- SX3221_rep3[,c(1,2,4,3)]

write.table(x = SX3221_rep3, 
            file = "SX3221_rep3_substracted.bdg", 
            row.names = F, 
            col.names = F, 
            quote = F, 
            sep = '\t')

