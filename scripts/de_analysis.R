#!/usr/bin/env Rscript

# load library
library("tidyverse")
library("Rsubread")
library("edgeR")
library("optparse")

#set ggplot preference
theme_set(theme_bw(base_size = 20)+
  theme(panel.grid = element_blank()))
colpal <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", 
            "#CC79A7", "#0072B2",  "#F0E442", "#999999")
    
options(ggplot2.discrete.colour=colpal, 
        ggplot2.discrete.fill=colpal,
        ggplot2.continuous.colour=colpal,
        ggplot2.continuous.fill=colpal)

#add options 
option_list <- list( 
    make_option(c("-b", "--basedir"), help="Base directory"),
    make_option(c("-n", "--annot"), help="Combined annotation"),
    make_option(c("-s", "--destab"), help="Design matrix of the analysis -> Anims, group"),
    make_option(c("-o", "--outfile"), help="Path to DE results"))

opt <- parse_args(OptionParser(option_list=option_list))
basedir  <- opt$basedir
gtfile  <- opt$annot
destab <- opt$design
outfile  <- opt$outfile


datcl  <- read.table(destab,header=FALSE,stringsAsFactors = FALSE)
colnames(datcl)  <- c("anims","treat")

bamlist  <- paste0(basedir,datanims$anims,".bam")

# Create count matrix
fc <- featureCounts(bamlist, annot.ext=gtfile,
                    isGTFAnnotationFile=TRUE, 
                    isPaired=FALSE,countMultiMappingReads=FALSE, nthreads=10)

# Create DGE object 

immun  <- as.factor(datcl$treat)
y  <- DGEList(rc,group=immun)


# Filter non-zero expression
y_nz <- y[rowSums(y$counts) > 0,]

# filter at least in 8 samples
y_fil <- y_nz[rowSums(cpm(y_nz)>=1) >=8,]

# replace DE with only filtered genes

y$counts  <- y_fil$counts

# calculate normalization factor
y <- calcNormFactors(y)

# create design matrix for comparison

design <- model.matrix(~immun)
y <- estimateDisp(y,design)
rownames(design)  <- colnames(y)


# Estimate dispersion matrix
y <- estimateDisp(y,design)

# Fit GLM on all filtered genes
fit <- glmQLFit(y, design)

# Test DE on all filtered genes
qlf <- glmQLFTest(fit)

# get the DE results
sortres <- topTags(object = qlf,n=nrow(qlf$table), adjust.method = "BH")

# DE results with Adj FDR <= 0.05
resfil  <- sortres$table  %>% filter(FDR <= 0.05)

# Write DE output results
write.table(resfil, file=outfile, sep="\t",row.names=FALSE, quote = FALSE)


