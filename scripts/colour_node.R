#!/usr/bin/env Rscript

library(tidyverse)

args  <- commandArgs(trailingOnly = TRUE)
datcore  <- read.table(paste0("remap/", args[1],"_coverage.tsv"), 
		  header=TRUE,
		  stringsAsFactors = FALSE)
assemb  <- args[2:length(args)]
colnames(datcore) <- c("nodeid","nodelen","chromo","pos","rrank",assemb)

##if has coverage put 1 in the matrix
datmat <- datcore[,colnames(datcore) %in% assemb] 
datmat[datmat<0] <- 0
datmat[datmat>0] <- 1


#add colour for nodes which already 
#defined from the rank
datlen <- read.table(paste0("graph/",args[1],"_graph_len.tsv"),
                     header=FALSE,stringsAsFactors = FALSE)

##get id match of the breed

datid <- data.frame(rankid=seq(0,length(assemb)-1),
                    breed=assemb,stringsAsFactors = FALSE)

colnames(datlen) <- c("nodeid","conlen","chromo","pos","rrank")

for (i in seq(1,nrow(datlen))){
  rankanims <- datlen[i,"rrank"]
  brid <- datid[datid$rankid==rankanims,"breed"]
  datmat[i,brid] <- 1
}


#extract colour from each node

datmat2 <- datmat 
datlab <- datmat 

colnode <- rep("0",nrow(datlab))
for (i in seq(1,nrow(datmat2))){
  labcol <- colnames(datmat2[i,which(datmat2[i,]==1),drop=FALSE]) %>% 
    paste0(collapse = ",")
  colnode[i] <- labcol
}

#output file
datout <- data.frame(nodeid=datcore$nodeid,colnode=colnode)

##node colour
write.table(datout,file=paste0("analysis/colour_node/",args[1],"_nodecol.tsv"),quote = FALSE,row.names = FALSE,col.names = TRUE)

#output as matrix
datmat$nodeid  <- datcore$nodeid
datmat  <- datmat  %>% select(nodeid, everything())
write.table(datmat,file=paste0("analysis/colour_node/",args[1],"_nodemat.tsv"),quote = FALSE,row.names = FALSE,col.names = TRUE)

