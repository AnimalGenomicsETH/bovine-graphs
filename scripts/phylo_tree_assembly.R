#!/usr/bin/env Rscript
### Given mash distance construct phylogenetic tree from it

#loading library 
library("tidyverse")
library("ape")

disfile  <- snakemake@input[[1]]
datdis  <- read.table(disfile,header=FALSE, stringsAsFactors =FALSE)


#rename the header
colnames(datdis)  <- c("anim1","anim2","distr","comp4","comp5")


# give correct assembly name

datdis$anim1c  <- str_extract(pattern ="([A-Z][A-Za-z]+)", string=datdis$anim1)
datdis$anim2c  <- str_extract(pattern ="([A-Z][A-Za-z]+)", string=datdis$anim2)

# make distance into a wide matrix 
datsel  <- datdis  %>% select(anim1c,anim2c, distr)
datwide  <- datsel  %>% pivot_wider(names_from = anim2c, values_from = distr)

datmat  <- as.matrix(datwide  %>% select(-anim1c))
rownames(datmat)  <- datwide$anim1c

outmat  <-  snakemake@output[[1]]

write.table(datmat,
	    file=outmat, sep="\t", quote = FALSE)

# select outgroup from the distance matrix
ref=snakemake@params[["ref"]]
outgroup  <- which.max(datmat[ref,])  %>% names()
print(outgroup)

# apply neighbor joining 

tr  <- nj(datmat)

# visualize the tree

outviz  <- snakemake@output[[2]]

pdf(outviz, width=12, height=10)

plot.phylo(root(tr, outgroup=outgroup), cex=2, edge.width=2)
axisPhylo(backward=FALSE, cex.axis=2)

dev.off()

