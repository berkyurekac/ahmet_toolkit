setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/RBP9/5dependent/piRNAtargets/bowtie/v_alignment")

a <- read.delim("piRNAs_new.txt", header = TRUE)
head(a)


b <- a[,c(1,7,9,10,12)]
colnames(b) <- c("genes", "shuffled_MM", "shuffled", "Sig_MM", "Sig")
head(b)

c <- subset(b, b$Sig_MM > b$Sig)

d <- subset(c, c$Sig_MM > c$shuffled_MM)

sig <- data.frame(d$genes,log2(d$Sig_MM/(d$shuffled_MM+1)))


shuffleMM <- read.delim("shuffled_MM2_pearl.sorted.HTSeqCount", header = FALSE, stringsAsFactors = FALSE)
shuffleNO <- read.delim("shuffled_NO.sorted.HTSeqCount", header = FALSE, stringsAsFactors = FALSE)
sigMM <- read.delim("sig_list_MM2_pearl.sorted.HTSeqCount", header = FALSE, stringsAsFactors = FALSE)
sigNO <- read.delim("sig_list_NO.sorted.HTSeqCount", header = FALSE, stringsAsFactors = FALSE)

HT <- data.frame(shuffleMM[,1:2], shuffleNO$V2, sigMM$V2, sigNO$V2)
colnames(HT) <- c("genes", "shuffled_MM", "shuffled", "Sig_MM", "Sig")
head(HT)

HTsig <- subset(HT, HT$Sig_MM > HT$Sig)

sig_list <- subset(sig, sig$log2.d.Sig_MM..d.shuffled_MM...1.. >= 1)
write.csv(sig_list, "rpb-9(mj261)_piRNAtargets.csv")


diff <- read.csv("lisa_ribozero_allgenes_padj.csv", header = TRUE, stringsAsFactors = FALSE)
diff2 <- read.csv("diffexpr-results.csv", header = TRUE, stringsAsFactors = FALSE)
overlap <- subset(diff, diff$X %in% c$genes)

overlapup <- subset(overlap, overlap$threshold == "A" )
overlapdown <- subset(overlap, overlap$threshold == "B" )












