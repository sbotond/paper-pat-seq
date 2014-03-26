#!/usr/bin/env Rscript

library(ggplot2)

# Deal with command line arguments directly (as optparse cannot handle multiple values):
args            <-list()
raw_args        <- commandArgs(trailingOnly = TRUE)
args$rep        <- raw_args[1]
args$in_files   <- raw_args[2:length(raw_args)]

pdf(args$rep)

lab<-list(
"runs_ma"       ="WT mean tail",
"runs_mb"       ="MUT mean tail",
"expr_gtail_a"  ="gtail coverage",
"expr_gtail_b"  ="gtail coverage"
)

# Function to threshold by coverage:
cutoff_plot <- function(d1, d2,c1,c2, limit, n1, n2) {
    i <- c1 > limit & c2 > limit & d1 > 0 & d2 > 0
    x <- d1[i]
    y <- d2[i]
    prop_cor <- cor(d1[i], d2[i],use="complete.obs")
    p <- ggplot(data.frame(x,y),aes(x=x,y=y)) + geom_point(size=1.1,alpha=0.4) + geom_smooth(method=lm)
    p<- p + ggtitle(sprintf("Tail lengths at coverage limit %d (r=%g)", l, prop_cor))
    p<- p + xlab(n1) + ylab(n2)
    print(p)
    return(prop_cor)
}

# Normalise rows and multiply with tail lengths:
fix_rows <- function(d){
    x       <- rep(c(0:(dim(d)[2]-1)),dim(d)[1])
    dim(x)  <- c(dim(d)[2],dim(d)[1])
    x       <- t(x)
    d       <- (d/rowSums(d)) * x
    return(d)
}

# Parse input file name:
parse_name <-function(s){
    x<-strsplit(s,"_|\\.")[[1]]
    return(x[length(x)-1])
}


limits  <-c(1,400,800)
in_files<-args$in_files

for(t1 in in_files) {
    for(t2 in in_files) {
    if(t1 == t2){ next }

    for(l in limits) {
       rd1  <- read.table(t1,header=TRUE,row.names=1)
       rd2  <- read.table(t2,header=TRUE,row.names=1)
       d1<-as.numeric(rowSums( fix_rows( rd1 ) ))
       d2<-as.numeric(rowSums( fix_rows( rd2 ) ))
       c1<-as.numeric(rowSums( rd1 ))
       c2<-as.numeric(rowSums( rd2 ))
       cutoff_plot(d1,d2,c1,c2, l, parse_name(t1), parse_name(t2)) 
    } 

    }
}

