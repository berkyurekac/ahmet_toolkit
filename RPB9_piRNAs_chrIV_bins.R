########
#this script was created by Ahmet Can Berkyurek on 18/04/2020, to analyze piRNA abundace along chrIV in wt, mutant and rescue!
########

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



######
# introduce raw coverage from mapped bam files, 5'dependent libraries
######

setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/RBP9/5dependent")

wt1 <- "isa_10_5dependentAligned.sortedByCoord.out_subset.bam"
## Coverage of a GAlignments object:
wt1_gal <- readGAlignments(wt1)
wt1_cvg <- coverage(wt1_gal)
wt1_cvg

wt3 <- "isa_12_5dependentAligned.sortedByCoord.out_subset.bam"
## Coverage of a GAlignments object:
wt3_gal <- readGAlignments(wt3)
wt3_cvg <- coverage(wt3_gal)
wt3_cvg


mut1 <- "isa_13_5dependentAligned.sortedByCoord.out_subset.bam"
## Coverage of a GAlignments object:
mut1_gal <- readGAlignments(mut1)
mut1_cvg <- coverage(mut1_gal)
mut1_cvg


mut2 <- "isa_14_5dependentAligned.sortedByCoord.out_subset.bam"
## Coverage of a GAlignments object:
mut2_gal <- readGAlignments(mut2)
mut2_cvg <- coverage(mut2_gal)
mut2_cvg


mut3 <- "isa_15_5dependentAligned.sortedByCoord.out_subset.bam"
## Coverage of a GAlignments object:
mut3_gal <- readGAlignments(mut3)
mut3_cvg <- coverage(mut3_gal)
mut3_cvg


rescue1 <- "isa_16_5dependentAligned.sortedByCoord.out_subset.bam"
## Coverage of a GAlignments object:
rescue1_gal <- readGAlignments(rescue1)
rescue1_cvg <- coverage(rescue1_gal)
rescue1_cvg


rescue2 <- "isa_17_5dependentAligned.sortedByCoord.out_subset.bam"
## Coverage of a GAlignments object:
rescue2_gal <- readGAlignments(rescue2)
rescue2_cvg <- coverage(rescue2_gal)
rescue2_cvg


rescue3 <- "isa_18_5dependentAligned.sortedByCoord.out_subset.bam"
## Coverage of a GAlignments object:
rescue3_gal <- readGAlignments(rescue3)
rescue3_cvg <- coverage(rescue3_gal)
rescue3_cvg

factor1 <- 2634941	
factor3 <- 1813824	
factor4 <- 9000499	
factor5 <- 8124378
factor6 <- 7396277	
factor7 <- 10924014	
factor8 <- 7661417
factor9 <- 3782503


#isolate ChrIV coverage for wt

vector1 <- as.vector(wt1_cvg[[4]])
vector3 <- as.vector(wt3_cvg[[4]])


a = as.data.frame(vector1)
c = as.data.frame(vector3)

# x axis
bps = 1:17493829

#y axis
coveragewt = data.frame(bps,a,c)

coveragewt$vector1 <- (coveragewt$vector1/factor1)*1e6
coveragewt$vector3 <- (coveragewt$vector3/factor3)*1e6



#isolate ChrIV coverage for mutant

vector4 <- as.vector(mut1_cvg$chrIV)
vector5 <- as.vector(mut2_cvg$chrIV)
vector6 <- as.vector(mut3_cvg$chrIV)


d = as.data.frame(vector4)
e = as.data.frame(vector5)
f = as.data.frame(vector6)

# x axis

#y axis
coveragemut = data.frame(bps,d,e,f)

coveragemut$vector4 <- (coveragemut$vector4/factor4)*1e6
coveragemut$vector5 <- (coveragemut$vector5/factor5)*1e6
coveragemut$vector6 <- (coveragemut$vector6/factor6)*1e6




#isolate ChrIV coverage for rescue

vector7 <- as.vector(rescue1_cvg$chrIV)
vector8 <- as.vector(rescue2_cvg$chrIV)
vector9 <- as.vector(rescue3_cvg$chrIV)


g = as.data.frame(vector7)
h = as.data.frame(vector8)
i = as.data.frame(vector9)

# x axis

#y axis
coveragerescue = data.frame(bps,g,h,i)

coveragerescue$vector7 <- (coveragerescue$vector7/factor4)*1e6
coveragerescue$vector8 <- (coveragerescue$vector8/factor5)*1e6
coveragerescue$vector9 <- (coveragerescue$vector9/factor6)*1e6



######
# compare wt/wt ,mut/wt, rescue/wt, do log2 transformation
######


merged <- data.frame(coveragewt, coveragemut, coveragerescue)

merged$meanwt <- rowMeans(merged[,2:3])
merged$meanmut <- rowMeans(merged[,5:7])
merged$meanrescue <- rowMeans(merged[,9:11])

head(merged)


merged_bar <- aggregate(merged[,1:18], by=list((seq(nrow(merged))-1) %/% 10000), FUN=mean)
values <- seq(1,17493829,10000)
merged_barframe <- data.frame(values,merged_bar)
head(merged_barframe)


x = ggplot(merged_barframe, aes(x=merged_barframe$values, y=log2(merged_barframe$meanrescue+1))) + geom_bar(stat="identity",color = "#0000FF") + scale_y_continuous(limits=c(-5,5)) + scale_x_continuous(limits=c(0,17493829)) +  guides(fill=FALSE) + xlab("ChrIV") + ylab('Reads per million') + theme_classic()

x2 = x + geom_bar(data=merged_barframe, aes(x=merged_barframe$values, y=-log2(merged_barframe$meanmut+1)), color="#FF0000", stat="identity")

x3 = x2 + geom_bar(data=merged_barframe, aes(x=merged_barframe$values, y=log2(merged_barframe$rescueoverwt)), color="blue", stat="identity")


x4 <- x3 + coord_flip() 
x5 <- x4 + scale_x_reverse()

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs")

ggsave("rescue_mut_RPM_piRNAs_rgb.pdf", x2, height = 5, width = 15, useDingbats=FALSE)

ggsave("mut_rescue_chrIV_new_coordflip2_all.pdf", x5, height = 15, width = 5, useDingbats=FALSE)




merged$wtoverwt <- ((merged$meanwt+1)/(merged$meanwt+1))

merged$mutoverwt <- ((merged$meanmut+1)/(merged$meanwt+1))

merged$rescueoverwt <- ((merged$meanrescue+1)/(merged$meanwt+1))

merged$wtovermut <- ((merged$meanwt+1)/(merged$meanmut+1))

head(merged)



merged_bar <- aggregate(merged[,15:18], by=list((seq(nrow(merged))-1) %/% 10000), FUN=mean)
values <- seq(1,17493829,10000)
merged_barframe <- data.frame(values,merged_bar)
head(merged_barframe)



x = ggplot(merged_barframe, aes(x=merged_barframe$values, y=log2(merged_barframe$mutoverwt))) + geom_bar(stat="identity",color = "Red") + scale_y_continuous(limits=c(-0.1,0.1)) + scale_x_continuous(limits=c(0,17493829)) +  guides(fill=FALSE) + xlab("ChrIV") + ylab('Reads per million') + theme_classic()

x2 = x + geom_bar(data=merged_barframe, aes(x=merged_barframe$values, y=log2(merged_barframe$wtovermut)), color="grey", stat="identity")

x3 = x2 + geom_bar(data=merged_barframe, aes(x=merged_barframe$values, y=log2(merged_barframe$rescueoverwt)), color="blue", stat="identity")


x4 <- x3 + coord_flip() 
x5 <- x4 + scale_x_reverse()

setwd("/Users/berkyurekahmetcan/Dropbox (ericmiskalab)/RPB-9_discussion/piRNAs")

ggsave("mut_rescue_chrIV_new_all.pdf", x3, height = 5, width = 15, useDingbats=FALSE)

ggsave("mut_rescue_chrIV_new_coordflip2_all.pdf", x5, height = 15, width = 5, useDingbats=FALSE)

############
# repeat the bining on the merged data frame
############


b <- c(rep("A",200), rep("B",200), rep("C",200), rep("D",200), rep("E",200), rep("F",200), rep("G",200),rep("H",200),rep("I",200),rep("J",200),rep("K",200),rep("L",200),rep("M",200),rep("N",200),rep("O",200),rep("P",200),rep("R",200),rep("S",200),rep("T",200),rep("U",201))

b <- as.vector(b)

merged$label <- b

head(merged)

x = ggplot(merged, aes(x=merged$label, y=log2(merged$mutoverwt))) + geom_boxplot(color = "Red") + scale_y_continuous(limits=c(-8,8)) +  guides(fill=FALSE) + xlab("ChrIV") + ylab('Reads per million') + theme_classic() + coord_flip()

x2 = x + geom_bar(data=merged, aes(x=merged$values, y=log2(merged$rescueoverwt)), color = "blue", stat = "identity") 

ggsave("mut_rescue_chrIV_new.pdf.pdf", x2, height = 10, width = 15, useDingbats=FALSE)




