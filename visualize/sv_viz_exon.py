#!/usr/bin/env python
import graph_viz_exon
from random import sample
from PyPDF2 import PdfFileMerger
from collections import defaultdict
import os
import re
import argparse


def main():
    graphtype = snakemake.wildcards["asb"]
    graphcomp = snakemake.params["assemb"]
    graph = graph_viz_exon.generate_edges(graphtype)
    nodeinf = graph_viz_exon.graph_info(graphtype)

    infile = open(snakemake.input["annot_file"])
    lines = infile.readlines()

    merger = PdfFileMerger(strict=False)

    # determine mutation types

    def get_sv(svfile, svtype="biallelic"):
        all_sv = defaultdict(list)
        with open(svfile) as file_:
            for line in file_:
                if svtype == "biallelic":
                    # svid, mutype, *_ = line.strip().split()
                    chromo, leftcoord, mutype, *_ = line.strip().split()
                    svid = f"{chromo}_{leftcoord}"
                if svtype == "multiallelic":
                    # 1_535561        1898    5138    AltIns  s35,s123493,s38
                    svid, startpos, stoppos, mutype, _ = line.strip().split()
                all_sv[svid].append(mutype)
        return all_sv

    all_sv = {**get_sv(snakemake.input["bialsv_file"]), **
              get_sv(snakemake.input["multisv_file"], svtype="multiallelic")}

    # b1_1_111515 1 111514 s1 s2 intergenic 1 1 111514 intergenic 1
    for line in lines[1:]:
        if re.search("exon|CDS", line):
            line_comp = line.strip().split()
            svid = line_comp[0]
            mutations = ",".join(all_sv[svid])
            left_coord = f"{line_comp[1]} {line_comp[2]}"
            start_node, stop_node = line_comp[3:5]
            feat_left = f" {line_comp[5]} {line_comp[6]}"
            right_coord = f"{line_comp[7]}_{line_comp[8]}"
            feat_right = f" {line_comp[9]} {line_comp[10]} "
            left_annot = left_coord + feat_left
            right_annot = right_coord + feat_right
            annot = f"{mutations}\n{left_annot}\n{right_annot}"
            graph_viz_exon.visualize_graph(graphtype, graphcomp, graph, nodeinf,
                                           start=start_node, stop=stop_node, annot=annot)
            merger.append(f"{graphtype}_{start_node}_{stop_node}.pdf")
            os.remove(f"{graphtype}_{start_node}_{stop_node}")
            os.remove(f"{graphtype}_{start_node}_{stop_node}.pdf")

    merger.write(snakemake.output[0])
    infile.close()


if __name__ == "__main__":
    main()
