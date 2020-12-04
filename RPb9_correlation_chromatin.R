# this script is written by Ahmet Can Berkyurek to check the correlation between total RNA-seq data of SET-32, SET-25 and RPB-9. 
# the question is whether or not RPB-9 dependent pathway correltes with chromatin pathways.

#introduce the necessary packages
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
library(Hmisc)
library(corrgram)

setwd("/Volumes/miska/Ahmet Berkyurek/sequencing/NGS/RNA-seq/alper/ribozero_total/differential_expression")

#first, introduce rpb-9 mutant results in data frame, including all genes. 

rpb9results <- read.csv("rpb9_vs_wt_ribozero_allgenes.csv", stringsAsFactors = FALSE)
rpb9results <- as.data.frame(rpb9results)
head(rpb9results)

rpb9resultsframe <- data.frame(rpb9results$X,rpb9results$log2FoldChange)
colnames(rpb9resultsframe) <- c("genes", "l2fc")
rownames(rpb9resultsframe) <- rpb9resultsframe$genes
head(rpb9resultsframe)

# introduce set-25 results in data frame, and reorder in the order of rpb-9 mutant results. 

set25results <- read.csv("set25_vs_wt_ribozero_allgenes.csv", stringsAsFactors = FALSE)
set25results <- as.data.frame(set25results)
head(set25results)

set25resultsframe <- data.frame(set25results$X,set25results$log2FoldChange)
colnames(set25resultsframe) <- c("genes","l2fc")
rownames(set25resultsframe) <- set25resultsframe$genes
head(set25resultsframe)

rpb9correlationset25 <- set25resultsframe[match(rpb9resultsframe$genes,set25resultsframe$genes),]
head(rpb9correlationset25)

# prepare the data frame for the correlation between set-25 and rpb-9 from the last matrix. then plot with ggplot!

set25correlation <- data.frame(rpb9resultsframe$l2fc,rpb9correlationset25$l2fc)
rownames(set25correlation) <- rownames(rpb9correlationset25)
head(set25correlation)

set25correlation$rpb9resultsframe.l2fc <- as.vector(as.numeric(set25correlation$rpb9resultsframe.l2fc))
set25correlation$rpb9correlationset25.l2fc <- as.vector(as.numeric(set25correlation$rpb9correlationset25.l2fc))
head(set25correlation)

set25cor <- cor(set25correlation$rpb9resultsframe.l2fc,set25correlation$rpb9correlationset25.l2fc, method = "spearman", use = "pairwise.complete.obs")
plot(set25correlation$rpb9resultsframe.l2fc,set25correlation$rpb9correlationset25.l2fc)

head(set25cor)


theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
                           color = "black", size = 1, 
                           xlab = "x", ylab = "y",
                           ticks = 10,
                           textsize = 3,
                           xlimit = max(abs(xvals),abs(yvals)),
                           ylimit = max(abs(yvals),abs(xvals)),
                           epsilon = max(xlimit,ylimit)/50){
  
  #INPUT:
  #xvals .- Values of x that will be plotted
  #yvals .- Values of y that will be plotted
  #xgeo  .- x intercept value for y axis
  #ygeo  .- y intercept value for x axis
  #color .- Default color for axis
  #size  .- Line size for axis
  #xlab  .- Label for x axis
  #ylab  .- Label for y axis
  #ticks .- Number of ticks to add to plot in each axis
  #textsize .- Size of text for ticks
  #xlimit .- Limit value for x axis 
  #ylimit .- Limit value for y axis
  #epsilon .- Parameter for small space
  
  
  #Create axis 
  xaxis <- data.frame(x_ax = c(-8, 8), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-8, 8))
  
  #Add axis
  theme.list <- 
    list(
      theme_void(), #Empty the current theme
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = xaxis),
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = yaxis),
      annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
      annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
      xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
      ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
    )
  
  #Add ticks programatically
  ticks_x <- round(seq(-8, 8, length.out = 4),2)
  ticks_y <- round(seq(-8, 8, length.out = 4),2)
  
  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:4){
    
    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
                        yt = c(xgeo + epsilon, xgeo - epsilon))
    
    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
                        yt = rep(ticks_y[k], 2))
    
    #Add ticks to geom line for x axis
    theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
                                             data = xtick, size = size, 
                                             color = color)
    
    #Add labels to the x-ticks
    theme.list[[nlist + 4*k-2]] <- annotate("text", 
                                            x = ticks_x[k], 
                                            y = ygeo - 2.5*epsilon,
                                            size = textsize,
                                            label = paste(ticks_x[k]))
    
    
    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
                                             data = ytick, size = size, 
                                             color = color)
    
    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
                                          x = xgeo - 2.5*epsilon, 
                                          y = ticks_y[k],
                                          size = textsize,
                                          label = paste(ticks_y[k]))
  }
  
  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)
}

breaks = seq(-8,8,2)


p <- ggplot(set25correlation) + theme_geometry(set25correlation$rpb9resultsframe.l2fc,set25correlation$rpb9correlationset25.l2fc) +
  geom_point(x=set25correlation$rpb9resultsframe.l2fc,y=set25correlation$rpb9correlationset25.l2fc, size=1, color = "black") + scale_x_continuous(limits = c(-8,8), breaks = breaks) + scale_y_continuous(limits = c(-8,8), breaks = breaks) + xlab("l2fc(set-25/wt)") + ylab("l2fc(rpb-9/wt)")  + theme_bw()


# introduce set-32 results in data frame, and reorder in the order of rpb-9 mutant results. 

set32results <- read.csv("set32_vs_wt_ribozero_allgenes.csv")
set32results <- as.data.frame(set32results)
head(set32results)

set32resultsframe <- data.frame(set32results$X,set32results$log2FoldChange)
colnames(set32resultsframe) <- c("genes","l2fc")
rownames(set32resultsframe) <- set32resultsframe$genes
head(set32resultsframe)

rpb9correlationset32 <- set32resultsframe[match(rpb9resultsframe$genes,set32resultsframe$genes),]
head(rpb9correlationset32)

# prepare the data frame for the correlation between set-25 and rpb-9 from the last matrix. then plot with ggplot!

set32correlation <- data.frame(rpb9resultsframe$l2fc,rpb9correlationset32$l2fc)
rownames(set32correlation) <- rownames(rpb9correlationset32)
head(set32correlation)

set32cor <- cor(set32correlation$rpb9resultsframe.l2fc,set32correlation$rpb9correlationset32.l2fc, method = "spearman", use = "pairwise.complete.obs")

head(set32cor)

theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
                           color = "black", size = 1, 
                           xlab = "x", ylab = "y",
                           ticks = 10,
                           textsize = 3,
                           xlimit = max(abs(xvals),abs(yvals)),
                           ylimit = max(abs(yvals),abs(xvals)),
                           epsilon = max(xlimit,ylimit)/50){
  
  #INPUT:
  #xvals .- Values of x that will be plotted
  #yvals .- Values of y that will be plotted
  #xgeo  .- x intercept value for y axis
  #ygeo  .- y intercept value for x axis
  #color .- Default color for axis
  #size  .- Line size for axis
  #xlab  .- Label for x axis
  #ylab  .- Label for y axis
  #ticks .- Number of ticks to add to plot in each axis
  #textsize .- Size of text for ticks
  #xlimit .- Limit value for x axis 
  #ylimit .- Limit value for y axis
  #epsilon .- Parameter for small space
  
  
  #Create axis 
  xaxis <- data.frame(x_ax = c(-8, 8), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-8, 8))
  
  #Add axis
  theme.list <- 
    list(
      theme_void(), #Empty the current theme
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = xaxis),
      geom_line(aes(x = x_ax, y = y_ax), color = color, size = size, data = yaxis),
      annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
      annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
      xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
      ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
    )
  
  #Add ticks programatically
  ticks_x <- round(seq(-8, 8, length.out = 4),2)
  ticks_y <- round(seq(-8, 8, length.out = 4),2)
  
  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:4){
    
    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
                        yt = c(xgeo + epsilon, xgeo - epsilon))
    
    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
                        yt = rep(ticks_y[k], 2))
    
    #Add ticks to geom line for x axis
    theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
                                             data = xtick, size = size, 
                                             color = color)
    
    #Add labels to the x-ticks
    theme.list[[nlist + 4*k-2]] <- annotate("text", 
                                            x = ticks_x[k], 
                                            y = ygeo - 2.5*epsilon,
                                            size = textsize,
                                            label = paste(ticks_x[k]))
    
    
    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
                                             data = ytick, size = size, 
                                             color = color)
    
    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
                                          x = xgeo - 2.5*epsilon, 
                                          y = ticks_y[k],
                                          size = textsize,
                                          label = paste(ticks_y[k]))
  }
  
  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)
}

breaks = seq(-8,8,2)


p <- ggplot(set32correlation) + theme_geometry(set32correlation$rpb9resultsframe.l2fc,set32correlation$rpb9correlationset32.l2fc) +
  geom_point(x=set32correlation$rpb9resultsframe.l2fc,y=set32correlation$rpb9correlationset32.l2fc, size=1, color = "black") + scale_x_continuous(limits = c(-8,8), breaks = breaks) + scale_y_continuous(limits = c(-8,8), breaks = breaks) + xlab("l2fc(set-32/wt)") + ylab("l2fc(rpb-9/wt)")  + theme_bw()








