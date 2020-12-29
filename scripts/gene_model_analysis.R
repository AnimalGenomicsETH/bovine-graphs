#!/usr/bin/env Rscript

library(tidyverse, warn.conflicts = FALSE)

infile  <-  snakemake@input[[1]]
outfile  <-  snakemake@output[[1]]

print(infile)
print(outfile)
datgff <- read.table(infile,
		     header=FALSE,
		     sep="\t",
		     stringsAsFactors = FALSE)

# name only the analyzed columns
colnames(datgff)[1:5] <- c("contig","progmod","seqfeat","startfeat","stopfeat")
colnames(datgff)[9] <- "geninfo"

#extract gene id
datgff$genid <- str_extract(pattern="g\\w+",string=datgff$geninfo)

# length of the sequence features
datgff$lenfeat <- with(datgff,abs(stopfeat-startfeat))

# filter gene with CDS less than 50 aa 

filgen <- datgff %>% 
	filter(seqfeat=="CDS") %>% 
	group_by(genid) %>% 
  	summarise(nc=sum(lenfeat)) %>% 
	filter(nc<150) %>% 
	select(genid) %>% 
	pull()

datgff2 <- datgff %>% filter(! genid %in% filgen) %>% droplevels()

seqfeatsum  <- data.frame(feature=character(), stats=character(), stringsAsFactors=FALSE)

##total predicted gene 
gentot  <- length(unique(datgff2$genid))

#from how many SVs 
svtot  <- length(unique(datgff2$contig))
statout  <- glue::glue('{gentot}({svtot})')

seqfeatsum  <- seqfeatsum  %>% add_row(feature="Total novel genes (distinct SV)", stats=statout)

# complete gene with TSS-TTS

tss <- datgff2  %>% filter(seqfeat=="tss") %>% select(genid) %>% pull() %>% unique() 
tts <- datgff2  %>% filter(seqfeat=="tts") %>% select(genid) %>% pull() %>% unique() 
combtt <- tss[tss %in% tts]


datgffc  <- datgff2  %>% filter(genid %in% combtt)

##total predicted gene 
gentot  <- length(unique(datgffc$genid))

#from how many SVs 
svtot  <- length(unique(datgffc$contig))

statout  <- glue::glue('{gentot}({svtot})')

seqfeatsum  <- seqfeatsum  %>% add_row(feature="Complete novel genes", stats=statout)

# further stats only consider complete genes

#' Add feature mean max and min to the sequence fatures catalogues
#' @param feature = feature name
#' @param feature = feature stat to summarize
#' @return updated features catalogues
addstat <- function (feature, statvec) {
	
	statout  <- glue::glue("{ round(mean(statvec), 2)} ({round(max(statvec), 2)} - {round(min(statvec), 2)})")

	seqfeatsum  <<- seqfeatsum  %>% add_row(feature=feature, stats=statout)
}

# transcript length
transcript  <- datgffc[datgffc$seqfeat=='transcript','lenfeat']
addstat(feature = "Transcript length", statvec = transcript)


# exon length
exon  <- datgffc[datgffc$seqfeat=="exon","lenfeat"]
addstat(feature = "Exon length", statvec = exon)


# exon length per gene
exongene  <- datgffc %>% filter(seqfeat=="exon") %>% group_by(genid) %>% 
  summarise(lexon=sum(lenfeat))  %>% pull(lexon) 

addstat(feature = "Exon length/gene", statvec=exongene)

# exon count per gene
exoncount <- datgffc %>% filter(seqfeat=="exon") %>% group_by(genid) %>% 
  summarise(cxon=n())  %>% pull(cxon)

addstat(feature = "Exon count/gene", statvec=exoncount)

# CDS stat
cds <- datgffc %>% filter(seqfeat=="CDS") %>% select(lenfeat) %>% pull()

addstat(feature = "CDS length", statvec=cds)

# CDS stat per gene
cdsgene <- datgffc %>% 
	filter(seqfeat=="CDS") %>% 
	group_by(genid) %>% 
	summarise(ncds=sum(lenfeat))  %>% pull(ncds)

addstat(feature = "CDS length/gene", statvec=cdsgene)

# Protein length is coding /3
addstat(feature = "protein length", statvec=cdsgene/3)

# save to outfile

write.table(seqfeatsum,
	    file=outfile,
	    sep=";",
           quote=FALSE,
	   row.names=FALSE)
