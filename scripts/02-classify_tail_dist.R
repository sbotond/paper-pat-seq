#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(mclust)

# Deal with command line arguments directly (as optparse cannot handle multiple values):
raw_args <- commandArgs(trailingOnly = TRUE)
NAME     <- raw_args[1]
MIN_COV  <- as.numeric(raw_args[2])
REP      <- raw_args[3]
INTABS   <- raw_args[4:length(raw_args)]
#

pdf(REP)

merge_slices <- function(pool){
    slices <-list()
    for(s in pool){
        d <-read.table(s,header=TRUE)  
        d$transcript<-NULL
        slices<-c(slices, list(d))
    }
    return(Reduce("+",slices))
}

cov_filter  <-function(d, min_cov){
    i <- rowSums(d) > min_cov
    return(d[i,])
}

norm_rows <- function(d){
    for(i in 1:nrow(d)){
        d[i,] <- d[i,]/sum(d[i,])     
    }
    return(d)
}

plot_densities <-function(data,title) {
    colnames(data)  <-NULL
    dt              <- t(data)    # transpose your data frame
    rownames(data)  <- NULL       # erase the automatic names
    dm              <- melt(dt)   # make data ggplot compatible
    p   <- ggplot(data=dm) + geom_smooth(se=FALSE,alpha=0.2,size=0.02, aes(x=Var1, y=value, colour=factor(Var2))) + theme(legend.position="none")
    p <- p + ggtitle(sprintf("%s %s", title, NAME)) + xlab("Length") + ylab("Proportion")
    print(p)
}

mclust_dists<-function(d) {
    m <- Mclust(d,G=1:35)
    return(m)
}

average_clusters <- function(d,m) {
    tmp <- list()
    for(i in seq(1,m$G)){
        cls <- d[m$classification == i,]
        tmp <- c(tmp,list(colMeans(cls)))
    }
    nd  <- tmp[[1]]
    for(j in seq(2,length(tmp))) {
        nd <- rbind(nd, tmp[[j]])
    }
    return(as.data.frame(nd,row.names=1:dim(nd)[1]))
}


td <- merge_slices(INTABS)
td <- norm_rows( cov_filter(td, MIN_COV))
plot_densities(td,title=sprintf("Tail run distributions (with coverage > %d) in", MIN_COV))
m  <- mclust_dists(td)
summary(m)
ac <- average_clusters(td,m)
plot_densities(ac,title=sprintf("Clustered tail run distributions (with coverage > %d) in", MIN_COV))

