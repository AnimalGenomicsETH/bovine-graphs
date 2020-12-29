#!/usr/bin/env Rscript

library("ggplot2", warn.conflicts = FALSE)
library("dplyr", warn.conflicts = FALSE)
library("tidyr", warn.conflicts = FALSE)
library("forcats", warn.conflicts = FALSE)
library("magrittr", warn.conflicts = FALSE)

# test comment
#' Given pangenome graph calculate the non-reference sequences grouped by assembly
#' @param nodemat: matrix of the node colour
#' @param graphlen: length of nodes in the graphs
#' @return data.frame of assembly and corresponding non-ref seq
calculate_nonref <- function(nodemat, graphlen, graphtype, outbase = ".") {
    datmat <- read.table(nodemat, header = TRUE, stringsAsFactors = FALSE)

    # add the nodelen
    datlen <- read.table(graphlen, header = FALSE, stringsAsFactors = FALSE)
    colnames(datlen) <- c("nodeid", "conlen", "chromo", "pos", "rrank")

    # join to get the node len in the matrix
    datmat$rrank <- datlen$rrank
    datmat$conlen  <- datlen$conlen

    # select only non-ref nodes
    # convert into long format
    datsel  <- datmat  %>% filter(rrank > 0)  %>% 
	    		   gather(key = "assembly", 
                                 value = "presence",
                                -rrank, -conlen, -nodeid) %>% 
			   filter(presence > 0)

    # get the non-ref count and length
    datout  <- datsel  %>% group_by(assembly)  %>% 
               		   summarize(no_segment=n(),
			   	     length_segment=sum(conlen))

   # do not include the reference assembly
    refassemb  <- colnames(datmat)[2]
    #datout   <- datout[-which.min(datout$no_segment),]
    datout   <- datout  %>% filter(assembly != refassemb)   


    outfile  <- file.path(outbase, paste0(graphtype, "_nonref_analysis.tsv"))

    write.table(datout, file=outfile, quote = FALSE, row.names=FALSE)

}

#' Analyze the non-reference sequences sharing among assemblies
#' @param: nodemat : matrix of node colour and graphlen:node len information 
#' @return: plot of the non-ref sharing count and length across assemblies
nonref_sharing_pattern <- function (nodemat, graphlen, graphtype, outbase=".") {
    datmat <- read.table(nodemat, header = TRUE, stringsAsFactors = FALSE)

    # add the nodelen
    datlen <- read.table(graphlen, header = FALSE, stringsAsFactors = FALSE)
    colnames(datlen) <- c("nodeid", "conlen", "chromo", "pos", "rrank")

    # join to get the node len in the matrix
    datmat$conlen  <- datlen$conlen
    # include only non-ref
    datsel  <- datmat[datlen$rrank>0,]

    # do not include reference column
    datsel  <- datsel[,-2]

    breeds  <- colnames(datsel  %>% select(-conlen, -nodeid))  %>% as.character()
    
    datsel$nonrefcol <- apply(datsel  %>% select(-conlen,-nodeid),
			      1,
	      		     function(x) paste(breeds[x >0], collapse=","))

    datsum  <- datsel  %>% group_by(nonrefcol)  %>% summarize(shared_count=n(),
    							      shared_len=sum(conlen))

    # plot 
    theme_set(theme_bw(base_size = 18)+
        theme(panel.grid = element_blank()))
    
    pl1 <- ggplot(datsum,aes(x=fct_reorder(nonrefcol,shared_len),y=shared_len/10^6))+
    geom_col(fill="#56B4E9",col="black",size=0.5)+
    geom_text(aes(label=round(shared_len/10^6,2)),hjust=-0.2,size=6)+
    scale_y_continuous(limits=c(0,max(datsum$shared_len)/10^6+5))+
    coord_flip()+
    labs(x="Bovine assemblies",
       y="Shared non-ref segment lengths (Mb)")


    pl2 <- ggplot(datsum,aes(x=fct_reorder(nonrefcol,shared_count),y=shared_count))+
    geom_col(fill="#009E73",col="black",size=0.5)+
    geom_text(aes(label=shared_count),hjust=-0.2,size=6)+
    scale_y_continuous(limits=c(0,max(datsum$shared_count)+5000))+
    coord_flip()+
    labs(x="Bovine assemblies",
       y="Shared non-ref segment counts")
    

    ggsave(pl1, filename = file.path(outbase, paste0(graphtype,"_nonref_shared_len.pdf")),width=12, height=10)
    ggsave(pl1, filename = file.path(outbase, paste0(graphtype,"_nonref_shared_len.png")),width=12, height=10)
    ggsave(pl2, filename = file.path(outbase, paste0(graphtype,"_nonref_shared_count.pdf")),width=12, height=10)
    ggsave(pl2, filename = file.path(outbase, paste0(graphtype,"_nonref_shared_count.png")),width=12, height=10)
    
}

