#!/usr/bin/env Rscript --silent

library(optparse)
library(ggplot2)
library(gridExtra)
library(plyr)
library(gtools)
library(RColorBrewer)

limits <- seq(0, 1000, 50)
x <- 1:length(limits)

row_mean <- function(x) { sum(0:(length(x)-1) * x)/sum(x) }

# Function to threshold by coverage
prop_cor <- function(d1, d2, prop, cut_prop, limit) {
  i <- d1[cut_prop] > limit & d2[cut_prop] > limit
  x <- d1[prop][i]
  y <- d2[prop][i]
  prop_cor <- cor(d1[prop][i], d2[prop][i],use="complete.obs")
  
  return(prop_cor)
}

rep_cor <- function(rep1, rep2, prop="mean", cut_prop="count") {
  vapply(limits, function(l) { prop_cor(rep1, rep2, prop=prop, cut_prop=cut_prop, limit=l) }, 0)
}

tech_cor <- function(files.df) {
  tech_reps = list()
  for (rep in files.df$file) {
    rep_raw = read.table(rep, header=TRUE, row.names=1)
    tech_reps[[rep]] = data.frame(mean=apply(rep_raw, 1, row_mean), count=apply(rep_raw, 1, sum))
  }
  
  perms <- permutations(n=length(tech_reps), r=2, v=names(tech_reps), repeats.allowed=F)
  
  cor_mat <- matrix(NA, nrow=0, ncol=length(limits))
  for (i in 1:nrow(perms)) { 
    perm <- perms[i, ]
    if (perm[1] > perm[2]) {
      next
    }
    cor_mat <- rbind(cor_mat, rep_cor(tech_reps[[perm[1]]], tech_reps[[perm[2]]]))
  }
  
  vals <- apply(cor_mat, 2, mean)
  
  data.frame(Dataset=files.df$set[1], Threshold=limits, Correlation=vals)
}


option_list <- list(
  make_option("--bio_replicates", default=NULL, help = "File containing of pairs of biological replicate datafiles."),
  make_option("--tech_replicates", default=NULL, help = "File containingsets of technical replicate datafiles."),
  make_option("--rep", default="replicate_correlations.pdf", help = "Report PDF.")
)
raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(OptionParser(option_list = option_list), args = raw_args)

# sets <- c(rep("MUT1", 4), rep("WT1", 4), rep("MUT2", 4), rep("WT2", 4))
# reps <- c("TEST_WT1_vs_MUT1_truns_MUT1A.tab", "TEST_WT1_vs_MUT1_truns_MUT1B.tab", "TEST_WT1_vs_MUT1_truns_MUT1C.tab", "TEST_WT1_vs_MUT1_truns_MUT1D.tab",
# "TEST_WT1_vs_MUT1_truns_WT1A.tab", "TEST_WT1_vs_MUT1_truns_WT1B.tab", "TEST_WT1_vs_MUT1_truns_WT1C.tab", "TEST_WT1_vs_MUT1_truns_WT1D.tab",
# "TEST_WT2_vs_MUT2_truns_MUT2A.tab", "TEST_WT2_vs_MUT2_truns_MUT2B.tab", "TEST_WT2_vs_MUT2_truns_MUT2C.tab", "TEST_WT2_vs_MUT2_truns_MUT2D.tab",
# "TEST_WT2_vs_MUT2_truns_WT2A.tab", "TEST_WT2_vs_MUT2_truns_WT2B.tab", "TEST_WT2_vs_MUT2_truns_WT2C.tab", "TEST_WT2_vs_MUT2_truns_WT2D.tab")
# tech_reps <- data.frame(Dataset=sets, file=reps)
# write.table(tech_reps, file="technical_replicates.tab", quote=F, row.names=F, sep="\t")
tech_reps <- read.table(args$tech_replicates, header=T)

all.tech.cors <- ddply(tech_reps, c("set"), tech_cor)
tech.cors.plot <- ggplot(aes(x=Threshold, y=Correlation, color=Dataset), data = all.tech.cors) + geom_line() + geom_point() + ylim(0, 1) +
  theme(legend.direction = "horizontal", legend.position = "bottom") + ggtitle("Mean correlation between technical replicates")

## 
## Biological replicates
##
bio_rep_design <- read.table(args$bio_replicates, header=T, as.is=T)

cor_mata <- matrix(NA, nrow=0, ncol=length(limits))
cor_matb <- matrix(NA, nrow=0, ncol=length(limits))
for (i in nrow(bio_rep_design)) {
  rep1  <- read.table(bio_rep_design$rep1[i], header=TRUE)
  rep2  <- read.table(bio_rep_design$rep2[i], header=TRUE)

  cor_mata <- rbind(cor_mata, rep_cor(rep1, rep2, prop="runs_ma", cut_prop="expr_gtail_a"))
  cor_matb <- rbind(cor_matb, rep_cor(rep1, rep2, prop="runs_mb", cut_prop="expr_gtail_b"))  
}

all.bio.cors <- data.frame(Dataset=c(rep("WT", length(limits)), rep("MUT", length(limits))), Threshold=rep(limits, 2),
                           Correlation=c(apply(cor_mata, 2, mean), apply(cor_matb, 2, mean)))
bio.cors.plot <- ggplot(aes(x=Threshold, y=Correlation, color=Dataset), data = all.bio.cors) + geom_line() + geom_point() + ylim(0, 1) +
  scale_color_discrete(breaks=c("MUT", "WT"), labels=c("MUT1 vs. MUT2", "WT1 vs. WT2")) +
  theme(legend.direction = "horizontal", legend.position = "bottom") + ggtitle("Correlation between biological replicates")

##
## Putting it together
##
ggsave(file=args$rep, plot=arrangeGrob(bio.cors.plot, tech.cors.plot, ncol=2), width=10, height=6, units="in")
