#!/usr/bin/env Rscript

library(optparse) # on CRAN
library(ggplot2)


# Deal with command line arguments:
option_list <- list(
   make_option("--d1",    default=NULL, help = "Data file from the first study."),
   make_option("--name1", default=NULL, help = "Name of the first dataset."),
   make_option("--tf1", default=NULL, help = "1: Transcript field."),
   make_option("--df1", default=NULL, help = "1: Data field."),
   make_option("--d2",    default=NULL, help = "Data file from the second study."),
   make_option("--name2", default=NULL, help = "Name of the second dataset."),
   make_option("--tf2", default=NULL, help = "2: Transcript field."),
   make_option("--df2", default=NULL, help = "2: Data field."),
   make_option("--rep", default="corr_with_sutdy.pdf", help = "Report PDF.")
   )
raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(OptionParser(option_list = option_list), args = raw_args)

pdf(args$rep)

raw1 <-read.table(args$d1,header=TRUE,stringsAsFactors=FALSE)
print(args$d2)
raw2 <-read.table(args$d2,header=TRUE,stringsAsFactors=FALSE)

lengths1 <-as.numeric(raw1[as.character(args$df1)][[1]])
lengths2 <-as.numeric(raw2[as.character(args$df2)][[1]])
trs1     <-as.character(raw1[as.character(args$tf1)][[1]])
trs2     <-as.character(raw2[as.character(args$tf2)][[1]])

names(lengths1) <- trs1
names(lengths2) <- trs2

tr <- intersect(trs1, trs2)

sc  <- cor(lengths1[tr], lengths2[tr],use="complete.obs")

p <- ggplot(data.frame(x=lengths1[tr],y=lengths2[tr]),aes(x=x,y=y)) + geom_point(size=1.1,alpha=0.4) + geom_smooth(method=lm)
p<- p + ggtitle(sprintf("Cross-study correlation (r=%g,n=%d)", sc, length(tr) ))
p<- p + xlab(args$name1) + ylab(args$name2)
print(p)

