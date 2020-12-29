#!/usr/bin/env Rscript

library(optparse)
snakemake@source("core_analysis.R")
snakemake@source("nonref_analysis.R")


#add options 
#option_list <- list( 
    #make_option(c("-g", "--graph"), help="Graph to calculate"))
#opt <- parse_args(OptionParser(option_list=option_list))
#graphtype  <- opt$graph
graphtype  <- snakemake@wildcards[["asb"]]


graphlen  <- paste0("graph/", graphtype, "_graph_len.tsv")
nodemat  <- paste0("analysis/colour_node/",graphtype, "_nodemat.tsv")
outbase  <- "analysis/core_nonref"

#run the core analysis
calculate_core(nodemat = nodemat, graphlen = graphlen, graphtype = graphtype, outbase = outbase)

pangenome_sampling(nodemat = nodemat, graphlen = graphlen, graphtype = graphtype, outbase = outbase)

#run the non-ref analysis

calculate_nonref(nodemat = nodemat, graphlen = graphlen, graphtype = graphtype, outbase = outbase)

nonref_sharing_pattern(nodemat = nodemat, graphlen = graphlen, graphtype = graphtype, outbase = outbase)
