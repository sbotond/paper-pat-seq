#!/usr/bin/env Rscript

library(optparse) # on CRAN


# Deal with command line arguments:
option_list <- list(
   make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
       help="Print extra output [default]"),
   make_option(c("-q", "--quietly"), action="store_false",
       dest="verbose", help="Print little output"),
   make_option(c("-c", "--count"), type="integer", default=5,
       help="Number of random normals to generate [default %default]",
       metavar="number"),
   make_option("--generator", default="rnorm",
       help = "Function to generate random deviates [default \"%default\"]"),
   make_option("--mean", default=0,
       help="Mean if generator == \"rnorm\" [default %default]"),
   make_option("--sd", default=1, metavar="standard deviation",
       help="Standard deviation if generator == \"rnorm\" [default %default]")
   )
raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(OptionParser(option_list = option_list), args = raw_args)

print(args)

