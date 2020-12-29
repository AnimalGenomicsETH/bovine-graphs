#!/usr/bin/env Rscript

library("dplyr", warn.conflicts = FALSE)
library("tidyr", warn.conflicts = FALSE)
library("ggplot2", warn.conflicts = FALSE)
library("magrittr", warn.conflicts = FALSE)


#' Given graph output the core and flexible genome size
#' @param graphlen: graphlen file in graph folder
#' @param nodemat: matrix of the node colour
#' @return the core (present in all assembly), flexible (not in all assembly)
#' flexible_shared at least in two assemblies, private: in a single assembly
calculate_core <- function(nodemat, graphlen, graphtype, outbase = ".") {
    datmat <- read.table(nodemat, header = TRUE, stringsAsFactors = FALSE)
    # add the number of colour in each node
    datmat$combres <- rowSums(datmat %>% select(-nodeid))

    # total number of assembly in the graph
    totassemb <- ncol(datmat) - 2

    # add the nodelen
    datlen <- read.table(graphlen, header = FALSE, stringsAsFactors = FALSE)
    colnames(datlen) <- c("nodeid", "conlen", "chromo", "pos", "rrank")

    # join to get the node len in the matrix
    datmat <- datmat %>% left_join(datlen, by = c("nodeid"))

    # calculate core genome length 
    core <- sum(datmat[datmat$combres == totassemb, "conlen"])
    flexible <- sum(datmat[datmat$combres < totassemb, "conlen"])
    flexible_shared <- sum(datmat[datmat$combres < totassemb & datmat$combres > 1, "conlen"])
    private <- sum(datmat[datmat$combres == 1, "conlen"])

    # add core and flexible node count 
    node_core <- length(datmat[datmat$combres == totassemb, "conlen"])
    node_flexible <- length(datmat[datmat$combres < totassemb, "conlen"])
    node_flexible_shared <- length(datmat[datmat$combres < totassemb & datmat$combres > 1, "conlen"])
    node_private <- length(datmat[datmat$combres == 1, "conlen"])

    outfile <- file.path(outbase, paste0(graphtype, "_core_analysis.tsv"))
    cat("core", "flexible", "flexible_shared", "private", "\n", file = outfile, append = TRUE)
    cat(core, flexible, flexible_shared, private, "\n", file = outfile, append = TRUE)
    cat(node_core, node_flexible, node_flexible_shared, node_private, file = outfile, append = TRUE)
}


#' Given pangenome graphs subsample and track the pangenome size changes
#' @params nodemat: matrix of the node colour
#' @params graphlen: info of node len
#' @return sampling results .tsv and plot as the sample increase in the graph
pangenome_sampling <- function(graphlen, nodemat, graphtype, outbase = ".") {
    ## experiment to get core genome as we add more genome into pangenome
    datexp <- read.table(nodemat, header = TRUE, stringsAsFactors = FALSE)


    # add the nodelen
    datlen <- read.table(graphlen, header = FALSE, stringsAsFactors = FALSE)
    colnames(datlen) <- c("nodeid", "conlen", "chromo", "pos", "rrank")

    # breeds in the graph
    breeds <- colnames(datexp %>% select(-nodeid))

    # number repetition in the sampling
    no_rep <- ncol(datexp) - 1

    # dataframe to store the sampling
    datpan <- data.frame(
        norep = numeric(),
        nosamp = numeric(),
        assemb = character(),
        core_gen = numeric(),
        flex_gen = numeric(),
        tot_gen = numeric()
    )

    datcon <- datlen %>% select(nodeid, conlen)

    # size of sampled
    for (nosamp in seq(1, length(breeds))) {
        # how many sampling repeated
        for (norep in seq(1, no_rep)) {
            selsamp <- as.character(sample(breeds, size = nosamp))
            # selected breeds
            selcol <- datexp[, colnames(datexp) %in% selsamp, drop = FALSE]
            # add contig length
            # it is ordered so we just add it from conlen
            # add number of colour in node
            selcol$comcol <- rowSums(selcol)
            # add contig len
            selcol$conlen <- datlen$conlen
            # core genome
            # core if shared in all of the member of population
            core_gen <- selcol[selcol$comcol == nosamp, "conlen"] %>% sum()
            # flex genome if shared less than the total of population
	    # not consider node not present 
	    if (nosamp == 1) 
		    flex_gen  <- 0 
	    else 
		    flex_gen <- selcol[selcol$comcol < nosamp & selcol$comcol > 0, "conlen"] %>% sum()
            # tot_genome if at least observed in a single breeds
            tot_gen <- selcol[selcol$comcol > 0, "conlen"] %>% sum()
            datpan <- rbind(datpan, data.frame(
                norep = norep,
                nosamp = nosamp,
                assemb = paste(selsamp, collapse = ","),
                core_gen = core_gen,
                flex_gen = flex_gen,
                tot_gen = tot_gen
            ))
        }
    }
    ## plot it
    panstat <- c(
        core_gen = "Core genome",
        flex_gen = "Flexible genome",
        tot_gen = "Total pangenome"
    )

    # convert to long format
    datlong <- datpan %>% gather(core_gen:tot_gen,
        key = "pan_stat", value = "val_stat"
    )

    # add mean line
    datmean <- datlong %>%
        group_by(nosamp, pan_stat) %>%
        summarise(mc = mean(val_stat))
    # main plotting command

    pl1 <- ggplot(datlong, aes(x = nosamp, y = val_stat / 10^6, group = nosamp)) +
        geom_jitter(width = 0.15) +
        geom_line(data = datmean, aes(y = mc / 10^6, group = 1), size = 1.5, col = "#009E73") +
        scale_x_continuous(breaks = seq(1, 6)) +
        facet_wrap(~pan_stat, nrow = 3, scale = "free_y", labeller = as_labeller(panstat)) +
        theme_bw(base_size = 18) +
        theme(
            strip.background = element_blank(),
            panel.grid = element_blank()
        ) +
        labs(
            x = "Number of genome(s) in the pangenome",
            y = "Size (Mb)"
        )

    # save plot
    ggsave(pl1,
        filename = file.path(outbase, paste0(graphtype, "_core_flex_sim.pdf")),
        width = 10, height = 12
    )

    ggsave(pl1,
        filename = file.path(outbase, paste0(graphtype, "_core_flex_sim.png")),
        width = 10, height = 12
    )


    # save the sampling result
    write.table(datpan,
        file = file.path(outbase, paste0(graphtype, "_core_flex_sim.tsv")),
        quote = FALSE, row.names = FALSE
    )
}
