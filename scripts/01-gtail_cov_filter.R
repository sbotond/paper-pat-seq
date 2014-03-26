#!/usr/bin/env Rscript

library(optparse) # on CRAN
library(ggplot2)


# Deal with command line arguments:
option_list <- list(
   make_option("--test1", default=NULL, help = "Transcript details from the first test."),
   make_option("--test1_out", default=NULL, help = "Filtered output for the first test."),
   make_option("--test2", default=NULL, help = "Transcript details from the second test."),
   make_option("--test2_out", default=NULL, help = "Filtered output for the second test."),
   make_option("--mc_wt", default=NULL, help = "Covergare threshold: WT."),
   make_option("--mc_mut", default=NULL, help = "Coverage threshold: MUT."),
   make_option("--rep", default="gtail_cov_filter.pdf", help = "Report PDF.")
   )
raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(OptionParser(option_list = option_list), args = raw_args)

pdf(args$rep)

t1  <- read.table(args$test1, header=TRUE)
t2  <- read.table(args$test2, header=TRUE)

lab<-list(
"runs_ma"       ="WT mean tail",
"runs_mb"       ="MUT mean tail",
"expr_gtail_a"  ="gtail coverage",
"expr_gtail_b"  ="gtail coverage"
)

# Function to threshold by coverage:
cutoff_plot <- function(d1, d2, prop, cut_prop, limit) {
    i <- d1[cut_prop] > limit & d2[cut_prop] > limit
    x <- d1[prop][i]
    y <- d2[prop][i]
    prop_cor <- cor(d1[prop][i], d2[prop][i],use="complete.obs")
    p <- ggplot(data.frame(x,y),aes(x=x,y=y)) + geom_point(size=1.1,alpha=0.4) + geom_smooth(method=lm)
    p<- p + ggtitle(sprintf("%s with a minimum %d %s (r=%g,nr=%d)",lab[prop],limit,lab[cut_prop],prop_cor,length(d1[prop][i])))
    p<- p + xlab(sprintf("TEST1: %s", lab[prop])) + ylab(sprintf("TEST2: %s", lab[prop]))
    print(p)
    return(prop_cor)
}

limits<-c(1,100,400,800,1000)

for(l in limits) {

cor_ma <- cutoff_plot(t1, t2, prop="runs_ma",cut_prop="expr_gtail_a", limit=l)
cor_mb <- cutoff_plot(t1, t2, prop="runs_mb",cut_prop="expr_gtail_b", limit=l)

}

# Report number of expressed transcripts:

# Filter by selected coverage thresholds:
mc_wt   <- as.numeric(args$mc_wt)
mc_mut  <- as.numeric(args$mc_mut)
cat(sprintf("Minimum G-tail coverage for WT %d\n", mc_wt))
cat(sprintf("Minimum G-tail coverage for MUT %d\n", mc_mut))

i1  <- t1$expr_gtail_a > mc_wt & t1$expr_gtail_b > mc_mut
i2  <- t2$expr_gtail_a > mc_wt & t2$expr_gtail_b > mc_mut
i   <- i1 & i2

t1_out <- t1[i, ]
t2_out <- t2[i, ]
cat( sprintf("Total number of transcripts after filtering: %d\n", length(which(i)) ))

# Write filtered output:
write.table(t1_out, file=args$test1_out,row.names=FALSE)
write.table(t2_out, file=args$test2_out,row.names=FALSE)

