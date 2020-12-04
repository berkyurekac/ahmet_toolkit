## this script is written by Ahmet Can Berkyurek on 31/10/2018 Gurdon Institute, Cambridge!
############## my purpose is to analyze RPB-9 RNA-seq data and create a graph which will support co-transcriptional silencing. 

library(GenomicAlignments)
library(Rsamtools)
library(ggplot2)
library(dplyr)
library(abind)
library(Rsamtools)
library(rtracklayer)

setwd("/Users/berkyurekahmetcan/Desktop/lisa")
###############v###wt elution libraries
elution_hrde1_1 <- "SLX-15940.NEBNext01Aligned.sortedByCoord.out.bam"
## Coverage of a GAlignments object:
elution_hrde1_1_gal <- readGAlignments(elution_hrde1_1)
elution_hrde1_1_cvg <- coverage(elution_hrde1_1_gal)
elution_hrde1_1_cvg

elution_hrde1_2 <- "SLX-15940.NEBNext02Aligned.sortedByCoord.out.bam"
## Coverage of a GAlignments object:
elution_hrde1_2_gal <- readGAlignments(elution_hrde1_2)
elution_hrde1_2_cvg <- coverage(elution_hrde1_2_gal)
elution_hrde1_2_cvg

elution_hrde1_3 <- "SLX-15940.NEBNext03Aligned.sortedByCoord.out.bam"
## Coverage of a GAlignments object:
elution_hrde1_3_gal <- readGAlignments(elution_hrde1_3)
elution_hrde1_3_cvg <- coverage(elution_hrde1_3_gal)
elution_hrde1_3_cvg


################wt input libraries
input_hrde1_1 <- "SLX-15940.NEBNext10Aligned.sortedByCoord.out.bam"
## Coverage of a GAlignments object:
input_hrde1_1_gal <- readGAlignments(input_hrde1_1)
input_hrde1_1_cvg <- coverage(input_hrde1_1_gal)
input_hrde1_1_cvg

input_hrde1_2 <- "SLX-15940.NEBNext11Aligned.sortedByCoord.out.bam"
## Coverage of a GAlignments object:
input_hrde1_2_gal <- readGAlignments(input_hrde1_2)
input_hrde1_2_cvg <- coverage(input_hrde1_2_gal)
input_hrde1_2_cvg

input_hrde1_3 <- "SLX-15940.NEBNext12Aligned.sortedByCoord.out.bam"
## Coverage of a GAlignments object:
input_hrde1_3_gal <- readGAlignments(input_hrde1_3)
input_hrde1_3_cvg <- coverage(input_hrde1_3_gal)
input_hrde1_3_cvg


countelutionhrde1_1 <- countBam(elution_hrde1_1)
countelutionhrde1_1value <- countelutionhrde1_1[[6]]

countelutionhrde1_2 <- countBam(elution_hrde1_2)
countelutionhrde1_2value <- countelutionhrde1_2[[6]]

countelutionhrde1_3 <- countBam(elution_hrde1_3)
countelutionhrde1_3value <- countelutionhrde1_3[[6]]

meanlibrarysizeelutionhrde1 <- c(countelutionhrde1_1value, countelutionhrde1_2value, countelutionhrde1_3value)
result.mean.elution.hrde1 <- mean(meanlibrarysizeelutionhrde1)
print(result.mean.elution.hrde1)


normalisationfactorhrde1_1 <- (countelutionhrde1_1value/result.mean.elution.hrde1)
normalisationfactorhrde1_2 <- (countelutionhrde1_2value/result.mean.elution.hrde1)
normalisationfactorhrde1_3 <- (countelutionhrde1_3value/result.mean.elution.hrde1)


#############library size normalisation factor for input, column 6, [[6]], gives the numbers of uniquely mapped reads
countinputhrde1_1 <- countBam(input_hrde1_1)
countinputhrde1_1value <- countinputhrde1_1[[6]]

countinputhrde1_2 <- countBam(input_hrde1_2)
countinputhrde1_2value <- countinputhrde1_2[[6]]

countinputhrde1_3 <- countBam(input_hrde1_3)
countinputhrde1_3value <- countinputhrde1_3[[6]]

meanlibrarysizeinputhrde1 <- c(countinputhrde1_1value, countinputhrde1_2value, countinputhrde1_3value)
result.mean.input.hrde1 <- mean(meanlibrarysizeinputhrde1)
print(result.mean.input.hrde1)


normalisationfactorinputhrde1_1 <- (countinputhrde1_1value/result.mean.input.hrde1)
normalisationfactorinputhrde1_2 <- (countinputhrde1_2value/result.mean.input.hrde1)
normalisationfactorinputhrde1_3 <- (countinputhrde1_3value/result.mean.input.hrde1)

#######################elution libraries
vector4 <- as.vector(elution_hrde1_1_cvg[[8]])
vector5 <- as.vector(elution_hrde1_2_cvg[[8]])
vector6 <- as.vector(elution_hrde1_3_cvg[[8]])

#[[8]] corresponds to the piRNAsensor!

g = as.data.frame(vector4[1:881])
h = as.data.frame(vector5[1:881])
i = as.data.frame(vector6[1:881])

# x axis
bps = 1:881

#y axis
coverage4 = data.frame(bps,g)
coverage5 = data.frame(bps,h)
coverage6 = data.frame(bps,i)

#take the average of every 10 rows in the coverage!
coverage_bar_hrde1 <- aggregate(coverage4[,2], by=list((seq(nrow(coverage4))-1) %/% 20), FUN=mean)
values <- seq(1,881,20)
coverage_bar_hrde1_frame <- data.frame(values,coverage_bar_hrde1$x)
p_bar_framehrde1 = ggplot(coverage_bar_hrde1_frame, aes(x=values, y=coverage_bar_hrde1.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate1", y = 'raw coverage') + theme_classic()

coverage2_bar_hrde1 <- aggregate(coverage5[,2], by=list((seq(nrow(coverage5))-1) %/% 20), FUN=mean)
values <- seq(1,881,20)
coverage2_bar_hrde1_frame <- data.frame(values,coverage2_bar_hrde1$x)
p_bar_frame2hrde1 = ggplot(coverage2_bar_hrde1_frame, aes(x=values, y=coverage2_bar_hrde1.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate2", y = 'raw coverage') + theme_classic()

coverage3_bar_hrde1 <- aggregate(coverage6[,2], by=list((seq(nrow(coverage6))-1) %/% 20), FUN=mean)
values <- seq(1,881,20)
coverage3_bar_hrde1_frame <- data.frame(values,coverage3_bar_hrde1$x)
p_bar_frame3hrde1 = ggplot(coverage3_bar_hrde1_frame, aes(x=values, y=coverage3_bar_hrde1.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate3", y = 'raw coverage') + theme_classic()

#raw coverage plots (optinal)
ggsave("raw-coveragewt_1.pdf", p_bar_framehrde1, height = 10, width = 15)
ggsave("raw-coveragewt_2.pdf", p_bar_frame2hrde1, height = 10, width = 15)
ggsave("raw-coveragewt_3.pdf", p_bar_frame3hrde1, height = 10, width = 15)


#divide coverage with library size factor for each
depth_norm_hrde1_1 <- coverage_bar_hrde1_frame$coverage_bar_hrde1.x/normalisationfactorhrde1_1
depth_norm_hrde1_2 <- coverage2_bar_hrde1_frame$coverage2_bar_hrde1.x/normalisationfactorhrde1_2
depth_norm_hrde1_3 <- coverage3_bar_hrde1_frame$coverage3_bar_hrde1.x/normalisationfactorhrde1_3

coverage_q4 = data.frame(values, depth_norm_hrde1_1)
coverage_q5 = data.frame(values,depth_norm_hrde1_2)
coverage_q6 = data.frame(values,depth_norm_hrde1_3)


#plot the normalized coverage plots for each replicate
q2_hrde1 = ggplot(coverage_q6, aes(x=values, y=depth_norm_hrde1_3)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1 (tm1200) replicate1", y = 'normalized coverage') + theme_classic()
q1_hrde1 = ggplot(coverage_q5, aes(x=values, y=depth_norm_hrde1_2)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1 (tm1200) replicate2", y = 'normalized coverage') + theme_classic()
q_hrde1 = ggplot(coverage_q4, aes(x=values, y=depth_norm_hrde1_1)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1 (tm1200) replicate3", y = 'normalized coverage') + theme_classic()

ggsave("coveragehrde1_1.pdf", q_hrde1, height = 10, width = 15)
ggsave("coveragehrde1_2.pdf", q1_hrde1, height = 10, width = 15)
ggsave("coveragehrde1_3.pdf", q2_hrde1, height = 10, width = 15)

####################### input libraries
vector4input <- as.vector(input_hrde1_1_cvg[[8]])
vector5input <- as.vector(input_hrde1_2_cvg[[8]])
vector6input <- as.vector(input_hrde1_3_cvg[[8]])

#[[8]] corresponds to the piRNAsensor!

j = as.data.frame(vector4input[1:881])
k = as.data.frame(vector5input[1:881])
l = as.data.frame(vector6input[1:881])

# x axis
bps = 1:881

#y axis
coverage4input = data.frame(bps,j)
coverage5input = data.frame(bps,k)
coverage6input = data.frame(bps,l)

#take the average of every 10 rows in the coverage!
coverage_bar_hrde1_input <- aggregate(coverage4input[,2], by=list((seq(nrow(coverage4input))-1) %/% 20), FUN=mean)
values <- seq(1,881,20)
coverage_bar_hrde1_input_frame <- data.frame(values,coverage_bar_hrde1_input$x)
p_bar_framehrde1_input = ggplot(coverage_bar_hrde1_input_frame, aes(x=values, y=coverage_bar_hrde1_input.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate1 input", y = 'raw coverage') + theme_classic()

coverage2_bar_hrde1_input <- aggregate(coverage5input[,2], by=list((seq(nrow(coverage5input))-1) %/% 20), FUN=mean)
values <- seq(1,881,20)
coverage2_bar_hrde1_input_frame <- data.frame(values,coverage2_bar_hrde1_input$x)
p_bar_frame2hrde1_input = ggplot(coverage2_bar_hrde1_input_frame, aes(x=values, y=coverage2_bar_hrde1_input.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate2 input", y = 'raw coverage') + theme_classic()

coverage3_bar_hrde1_input <- aggregate(coverage6input[,2], by=list((seq(nrow(coverage6input))-1) %/% 20), FUN=mean)
values <- seq(1,881,20)
coverage3_bar_hrde1_input_frame <- data.frame(values,coverage3_bar_hrde1_input$x)
p_bar_frame3hrde1_input = ggplot(coverage3_bar_hrde1_input_frame, aes(x=values, y=coverage3_bar_hrde1_input.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate3 input", y = 'raw coverage') + theme_classic()

#raw coverage plots (optinal)
ggsave("raw-coveragehrde1input_1.pdf", p_bar_framehrde1_input, height = 10, width = 15)
ggsave("raw-coveragehrde1input_2.pdf", p_bar_frame2hrde1_input, height = 10, width = 15)
ggsave("raw-coveragehrde1input_3.pdf", p_bar_frame3hrde1_input, height = 10, width = 15)

#divide coverage with library size factor for each
depth_norm_inputhrde1_1 <- coverage_bar_hrde1_input_frame$coverage_bar_hrde1_input.x/normalisationfactorinputhrde1_1
depth_norm_inputhrde1_2 <- coverage2_bar_hrde1_input_frame$coverage2_bar_hrde1_input.x/normalisationfactorinputhrde1_2
depth_norm_inputhrde1_3 <- coverage3_bar_hrde1_input_frame$coverage3_bar_hrde1_input.x/normalisationfactorinputhrde1_3

coverage_q4input = data.frame(values, depth_norm_inputhrde1_1)
coverage_q5input = data.frame(values,depth_norm_inputhrde1_2)
coverage_q6input = data.frame(values,depth_norm_inputhrde1_3)


#plot the normalized coverage plots for each replicate
q6input = ggplot(coverage_q6input, aes(x=values, y=depth_norm_inputhrde1_3)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate3 input", y = 'normalized coverage') + theme_classic()
q5input = ggplot(coverage_q5input, aes(x=values, y=depth_norm_inputhrde1_2)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate2 input", y = 'normalized coverage') + theme_classic()
q4input = ggplot(coverage_q4input, aes(x=values, y=depth_norm_inputhrde1_1)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate1 input", y = 'normalized coverage') + theme_classic()

ggsave("coveragehrde1_1inputnormalized.pdf", q4input, height = 10, width = 15)
ggsave("coveragehrde1_2inputnormalized.pdf", q5input, height = 10, width = 15)
ggsave("coveragehrde1_3inputnormalized.pdf", q6input, height = 10, width = 15)

## fold change between elution and input hrde1

bps = 1:881

depth_foldchange_hrde1_1 = coverage_q4$depth_norm_hrde1_1/(coverage_q4input$depth_norm_inputhrde1_1 + 0.1)
coverage_foldchange_hrde1_1 <- data.frame(values,depth_foldchange_hrde1_1)
colnames(coverage_foldchange_hrde1_1) <- c('values', 'hrde1_1_norm_cov')


depth_foldchange_hrde1_2 = coverage_q5$depth_norm_hrde1_2/(coverage_q5input$depth_norm_inputhrde1_2 + 0.1)
coverage_foldchange_hrde1_2 <- data.frame(values,depth_foldchange_hrde1_2)
colnames(coverage_foldchange_hrde1_2) <- c('values', 'hrde1_2_norm_cov')


depth_foldchange_hrde1_3 = (coverage_q6$depth_norm_hrde1_3/(coverage_q6input$depth_norm_inputhrde1_3 + 0.1))
coverage_foldchange_hrde1_3 <- data.frame(values,depth_foldchange_hrde1_3)
colnames(coverage_foldchange_hrde1_3) <- c('values', 'hrde1_3_norm_cov')


fchrde_1_1 = ggplot(coverage_foldchange_hrde1_1, aes(x=values, y=hrde1_1_norm_cov)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881, 1281, 1401, 1661)) + scale_y_continuous(limits=c(0,120)) + labs(title = "hrde-1(tm1200) replicate-1", y = 'Fold Change (Elution/Input)') + theme_classic()
fchrde_1_2 = ggplot(coverage_foldchange_hrde1_2, aes(x=values, y=hrde1_2_norm_cov)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881, 1281, 1401, 1661)) + scale_y_continuous(limits=c(0,120)) + labs(title = "hrde-1(tm1200) replicate-2", y = 'Fold Change (Elution/Input)') + theme_classic()
fchrde_1_3 = ggplot(coverage_foldchange_hrde1_3, aes(x=values, y=hrde1_3_norm_cov)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881, 1281, 1401, 1661)) + scale_y_continuous(limits=c(0,120)) + labs(title = "hrde-1(tm1200) replicate-3", y = 'Fold Change (Elution/Input)') + theme_classic()

ggsave("fold_hrde1_1.pdf",fchrde_1_1)
ggsave("fold_hrde1_2.pdf",fchrde_1_2)
ggsave("fold_hrde1_3.pdf",fchrde_1_3)


## normalisation to control! first combine all replicates and the the mean!
# install.packages('abind')

temp_arrayhrde1 <- abind(coverage_foldchange_hrde1_1, coverage_foldchange_hrde1_2, coverage_foldchange_hrde1_3, along=3)
coveragehrde1_combined <- apply(temp_arrayhrde1, 1:2, mean)
coveragehrde1_combined_frame <- as.data.frame(coveragehrde1_combined)
colnames(coveragehrde1_combined_frame) <- c('values', 'meancoverage4')
head(coveragehrde1_combined_frame)

fchrde1_combined = ggplot(coveragehrde1_combined_frame, aes(x=values, y=meancoverage4)) + geom_smooth(method = "loess") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,2)) + labs(title = "hrde-1(tm1200) replicates mean", y = 'Fold Change (Elution/Input)') + theme_classic()
ggsave("wtvsmut_combined.pdf", fchrde1_combined, height = 10, width = 15)



# repeat the same for a control gene such as ama-1
vector4 <- as.vector(elution_hrde1_1_cvg[[4]])
vector5 <- as.vector(elution_hrde1_2_cvg[[4]])
vector6 <- as.vector(elution_hrde1_3_cvg[[4]])

#[[8]] corresponds to the piRNAsensor!

g = as.data.frame(vector4[4248020:4258173])
h = as.data.frame(vector5[4248020:4258173])
i = as.data.frame(vector6[4248020:4258173])

# x axis
bps = 1:10154

#y axis
coverage4 = data.frame(bps,g)
coverage5 = data.frame(bps,h)
coverage6 = data.frame(bps,i)

#take the average of every 10 rows in the coverage!
coverage_bar_hrde1 <- aggregate(coverage4[,2], by=list((seq(nrow(coverage4))-1) %/% 20), FUN=mean)
values <- seq(1,10154,20)
coverage_bar_hrde1_frame <- data.frame(values,coverage_bar_hrde1$x)
p_bar_framehrde1 = ggplot(coverage_bar_hrde1_frame, aes(x=values, y=coverage_bar_hrde1.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate1", y = 'raw coverage') + theme_classic()

coverage2_bar_hrde1 <- aggregate(coverage5[,2], by=list((seq(nrow(coverage5))-1) %/% 20), FUN=mean)
values <- seq(1,10154,20)
coverage2_bar_hrde1_frame <- data.frame(values,coverage2_bar_hrde1$x)
p_bar_frame2hrde1 = ggplot(coverage2_bar_hrde1_frame, aes(x=values, y=coverage2_bar_hrde1.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate2", y = 'raw coverage') + theme_classic()

coverage3_bar_hrde1 <- aggregate(coverage6[,2], by=list((seq(nrow(coverage6))-1) %/% 20), FUN=mean)
values <- seq(1,10154,20)
coverage3_bar_hrde1_frame <- data.frame(values,coverage3_bar_hrde1$x)
p_bar_frame3hrde1 = ggplot(coverage3_bar_hrde1_frame, aes(x=values, y=coverage3_bar_hrde1.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate3", y = 'raw coverage') + theme_classic()

#raw coverage plots (optinal)
ggsave("raw-coveragewt_ama1.pdf", p_bar_framehrde1, height = 10, width = 15)
ggsave("raw-coveragewt_ama2.pdf", p_bar_frame2hrde1, height = 10, width = 15)
ggsave("raw-coveragewt_ama3.pdf", p_bar_frame3hrde1, height = 10, width = 15)


#divide coverage with library size factor for each
depth_norm_hrde1_1 <- coverage_bar_hrde1_frame$coverage_bar_hrde1.x/normalisationfactorhrde1_1
depth_norm_hrde1_2 <- coverage2_bar_hrde1_frame$coverage2_bar_hrde1.x/normalisationfactorhrde1_2
depth_norm_hrde1_3 <- coverage3_bar_hrde1_frame$coverage3_bar_hrde1.x/normalisationfactorhrde1_3

coverage_q4 = data.frame(values, depth_norm_hrde1_1)
coverage_q5 = data.frame(values,depth_norm_hrde1_2)
coverage_q6 = data.frame(values,depth_norm_hrde1_3)


#plot the normalized coverage plots for each replicate
q2_hrde1 = ggplot(coverage_q6, aes(x=values, y=depth_norm_hrde1_3)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1 (tm1200) replicate1", y = 'normalized coverage') + theme_classic()
q1_hrde1 = ggplot(coverage_q5, aes(x=values, y=depth_norm_hrde1_2)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1 (tm1200) replicate2", y = 'normalized coverage') + theme_classic()
q_hrde1 = ggplot(coverage_q4, aes(x=values, y=depth_norm_hrde1_1)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1 (tm1200) replicate3", y = 'normalized coverage') + theme_classic()

ggsave("coveragehrde1_ama1.pdf", q_hrde1, height = 10, width = 15)
ggsave("coveragehrde1_ama2.pdf", q1_hrde1, height = 10, width = 15)
ggsave("coveragehrde1_ama3.pdf", q2_hrde1, height = 10, width = 15)




####################### input libraries
vector4input <- as.vector(input_hrde1_1_cvg[[4]])
vector5input <- as.vector(input_hrde1_2_cvg[[4]])
vector6input <- as.vector(input_hrde1_3_cvg[[4]])

#[[8]] corresponds to the piRNAsensor!

j = as.data.frame(vector4input[4248020:4258173])
k = as.data.frame(vector5input[4248020:4258173])
l = as.data.frame(vector6input[4248020:4258173])

# x axis
bps = 1:10154

#y axis
coverage4input = data.frame(bps,j)
coverage5input = data.frame(bps,k)
coverage6input = data.frame(bps,l)

#take the average of every 10 rows in the coverage!
coverage_bar_hrde1_input <- aggregate(coverage4input[,2], by=list((seq(nrow(coverage4input))-1) %/% 20), FUN=mean)
values <- seq(1,10154,20)
coverage_bar_hrde1_input_frame <- data.frame(values,coverage_bar_hrde1_input$x)
p_bar_framehrde1_input = ggplot(coverage_bar_hrde1_input_frame, aes(x=values, y=coverage_bar_hrde1_input.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate1 input", y = 'raw coverage') + theme_classic()

coverage2_bar_hrde1_input <- aggregate(coverage5input[,2], by=list((seq(nrow(coverage5input))-1) %/% 20), FUN=mean)
values <- seq(1,10154,20)
coverage2_bar_hrde1_input_frame <- data.frame(values,coverage2_bar_hrde1_input$x)
p_bar_frame2hrde1_input = ggplot(coverage2_bar_hrde1_input_frame, aes(x=values, y=coverage2_bar_hrde1_input.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE)  + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate2 input", y = 'raw coverage') + theme_classic()

coverage3_bar_hrde1_input <- aggregate(coverage6input[,2], by=list((seq(nrow(coverage6input))-1) %/% 20), FUN=mean)
values <- seq(1,10154,20)
coverage3_bar_hrde1_input_frame <- data.frame(values,coverage3_bar_hrde1_input$x)
p_bar_frame3hrde1_input = ggplot(coverage3_bar_hrde1_input_frame, aes(x=values, y=coverage3_bar_hrde1_input.x), geom="histogram") + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE)  + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate3 input", y = 'raw coverage') + theme_classic()

#raw coverage plots (optinal)
ggsave("raw-coveragehrde1ama1_1.pdf", p_bar_framehrde1_input, height = 10, width = 15)
ggsave("raw-coveragehrde1ama1_2.pdf", p_bar_frame2hrde1_input, height = 10, width = 15)
ggsave("raw-coveragehrde1ama1_3.pdf", p_bar_frame3hrde1_input, height = 10, width = 15)

#divide coverage with library size factor for each
depth_norm_inputhrde1_1 <- coverage_bar_hrde1_input_frame$coverage_bar_hrde1_input.x/normalisationfactorinputhrde1_1
depth_norm_inputhrde1_2 <- coverage2_bar_hrde1_input_frame$coverage2_bar_hrde1_input.x/normalisationfactorinputhrde1_2
depth_norm_inputhrde1_3 <- coverage3_bar_hrde1_input_frame$coverage3_bar_hrde1_input.x/normalisationfactorinputhrde1_3

coverage_q4input = data.frame(values, depth_norm_inputhrde1_1)
coverage_q5input = data.frame(values,depth_norm_inputhrde1_2)
coverage_q6input = data.frame(values,depth_norm_inputhrde1_3)


#plot the normalized coverage plots for each replicate
q6input = ggplot(coverage_q6input, aes(x=values, y=depth_norm_inputhrde1_3)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate3 input", y = 'normalized coverage') + theme_classic()
q5input = ggplot(coverage_q5input, aes(x=values, y=depth_norm_inputhrde1_2)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE)  + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate2 input", y = 'normalized coverage') + theme_classic()
q4input = ggplot(coverage_q4input, aes(x=values, y=depth_norm_inputhrde1_1)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_y_continuous(limits=c(0,200)) + labs(title = "hrde-1(tm1200) replicate1 input", y = 'normalized coverage') + theme_classic()

ggsave("coverageama1_1inputnormalized.pdf", q4input, height = 10, width = 15)
ggsave("coverageama1_2inputnormalized.pdf", q5input, height = 10, width = 15)
ggsave("coverageama1_3inputnormalized.pdf", q6input, height = 10, width = 15)



## fold change between elution and input hrde1

bps = 1:10154

depth_foldchange_hrde1_1 = coverage_q4$depth_norm_hrde1_1/(coverage_q4input$depth_norm_inputhrde1_1 + 1)
coverage_foldchange_hrde1_1 <- data.frame(values,depth_foldchange_hrde1_1)
colnames(coverage_foldchange_hrde1_1) <- c('values', 'hrde1_1_norm_cov')


depth_foldchange_hrde1_2 = coverage_q5$depth_norm_hrde1_2/(coverage_q5input$depth_norm_inputhrde1_2 + 1)
coverage_foldchange_hrde1_2 <- data.frame(values,depth_foldchange_hrde1_2)
colnames(coverage_foldchange_hrde1_2) <- c('values', 'hrde1_2_norm_cov')


depth_foldchange_hrde1_3 = (coverage_q6$depth_norm_hrde1_3/(coverage_q6input$depth_norm_inputhrde1_3 + 1))
coverage_foldchange_hrde1_3 <- data.frame(values,depth_foldchange_hrde1_3)
colnames(coverage_foldchange_hrde1_3) <- c('values', 'hrde1_3_norm_cov')


fchrde_1_1 = ggplot(coverage_foldchange_hrde1_1, aes(x=values, y=hrde1_1_norm_cov)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881, 1281, 1401, 1661)) + scale_y_continuous(limits=c(0,120)) + labs(title = "hrde-1(tm1200) replicate-1", y = 'Fold Change (Elution/Input)') + theme_classic()
fchrde_1_2 = ggplot(coverage_foldchange_hrde1_2, aes(x=values, y=hrde1_2_norm_cov)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881, 1281, 1401, 1661)) + scale_y_continuous(limits=c(0,120)) + labs(title = "hrde-1(tm1200) replicate-2", y = 'Fold Change (Elution/Input)') + theme_classic()
fchrde_1_3 = ggplot(coverage_foldchange_hrde1_3, aes(x=values, y=hrde1_3_norm_cov)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881, 1281, 1401, 1661)) + scale_y_continuous(limits=c(0,120)) + labs(title = "hrde-1(tm1200) replicate-3", y = 'Fold Change (Elution/Input)') + theme_classic()

ggsave("fold_hrde1_1.pdf",fchrde_1_1)
ggsave("fold_hrde1_2.pdf",fchrde_1_2)
ggsave("fold_hrde1_3.pdf",fchrde_1_3)




## normalisation to control! first combine all replicates and the the mean!
# install.packages('abind')

temp_arrayhrde1 <- abind(coverage_foldchange_hrde1_1, coverage_foldchange_hrde1_2, coverage_foldchange_hrde1_3, along=3)
coveragehrde1_combined <- apply(temp_arrayhrde1, 1:2, mean)
coveragehrde1_combined_frame <- as.data.frame(coveragehrde1_combined)
colnames(coveragehrde1_combined_frame) <- c('values', 'meancoverage4')
head(coveragehrde1_combined_frame)

fchrde1_combined = ggplot(coveragehrde1_combined_frame, aes(x=values, y=meancoverage4)) + geom_smooth(method = "loess") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881)) + scale_y_continuous(limits=c(0,2)) + labs(title = "hrde-1(tm1200) replicates mean", y = 'Fold Change (Elution/Input)') + theme_classic()
ggsave("ama_1_wtvsmut_combined.pdf", fchrde1_combined, height = 10, width = 15)


## combine wt and input normalized covrage !

# elution
temp_array_hrde1normalized <- abind(coverage_q4, coverage_q5, coverage_q6, along=3)
coveragehrde1normalized_combined <- apply(temp_array_hrde1normalized, 1:2, mean)
coveragehrde1_combinednormalized_frame <- as.data.frame(coveragehrde1normalized_combined)
colnames(coveragehrde1_combinednormalized_frame) <- c('values', 'meancoverage5')
head(coveragehrde1_combinednormalized_frame)

fchrde1_combinednormalized = ggplot(coveragehrde1_combinednormalized_frame, aes(x=values, y=meancoverage5)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881, 1281, 1401, 1661)) + scale_y_continuous(limits=c(0,120)) + labs(title = "hrde-1(tm1200) replicates mean", y = 'normalized coverage') + theme_classic()
ggsave("hrde1_combinednormalized.pdf", fchrde1_combinednormalized)

# input
temp_array_hrde1normalizedinput <- abind(coverage_q4input, coverage_q5input, coverage_q6input, along=3)
coveragehrde1normalizedinput_combined <- apply(temp_array_hrde1normalizedinput, 1:2, mean)
coveragehrde1_combinednormalizedinput_frame <- as.data.frame(coveragehrde1normalizedinput_combined)
colnames(coveragehrde1_combinednormalizedinput_frame) <- c('values', 'meancoverage6')
head(coveragehrde1_combinednormalizedinput_frame)

fchrde1_combinednormalizedinput = ggplot(coveragehrde1_combinednormalizedinput_frame, aes(x=values, y=meancoverage6)) + geom_bar(colour="grey", stat="identity") +  guides(fill=FALSE) + scale_x_continuous(breaks = c(1, 881, 1281, 1401, 1661)) + scale_y_continuous(limits=c(0,120)) + labs(title = "hrde-1(tm1200) replicates mean", y = 'normalized coverage') + theme_classic()
ggsave("hrde1_combinednormalizedinput.pdf", fchrde1_combinednormalizedinput)

#  Now, let's try to do this for the 58 gene list. filter the bam file and try to scale these guys! 
# first, using bedtools coverage function, determine the coverage on these 58 genes. 
# introduce each bed file and then, combined them! this is going to be a bit challenging. I don't know how to aggregate them!

wt1 <- read.delim("SLX-15940.NEBNext01Aligned.sortedByCoord.out.bam.bed", header = FALSE, stringsAsFactors = FALSE)
wt2 <- read.delim("SLX-15940.NEBNext02Aligned.sortedByCoord.out.bam.bed", header = FALSE, stringsAsFactors = FALSE)
wt3 <- read.delim("SLX-15940.NEBNext03Aligned.sortedByCoord.out.bam.bed", header = FALSE, stringsAsFactors = FALSE)

wt1$V4 <- as.numeric(as.vector(wt1$V4))
wt2$V4 <- as.numeric(as.vector(wt2$V4))
wt3$V4 <- as.numeric(as.vector(wt3$V4))

wt1set <- (subset(wt1, wt1$V4 == "gpa-16"))
wt2set <- (subset(wt2, wt2$V4 == "gpa-16"))
wt3set <- (subset(wt3, wt3$V4 == "gpa-16"))

gpa <- data.frame(wt1set$V7,(wt1set$V4 + wt2set$V4 + wt3set$V4)/3)

head(wt1set)
tail(wt1set)

mut1 <- read.delim("SLX-15940.NEBNext10Aligned.sortedByCoord.out.bam.bed", header = FALSE)
mut2 <- read.delim("SLX-15940.NEBNext11Aligned.sortedByCoord.out.bam.bed", header = FALSE)
mut3 <- read.delim("SLX-15940.NEBNext12Aligned.sortedByCoord.out.bam.bed", header = FALSE)







