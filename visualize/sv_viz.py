#!/usr/bin/env python
import graph_viz
from random import sample
from PyPDF2 import PdfFileMerger
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--graphtype", help="graph to process")
    parser.add_argument("-c", "--graphcomp", help="assembly in the graph", nargs="+")
    parser.add_argument("-m", "--mutype", help="mutation type", type=str)
    return parser.parse_args()


def main():
    args = parse_args()
    graphtype = args.graphtype
    graphcomp = args.graphcomp
    mutype = args.mutype
    graph = graph_viz.generate_edges(graphtype)
    nodeinf = graph_viz.graph_info(graphtype)

    print(f"analysis/bubble/{graphtype}_{mutype}_sv.tsv")
    infile = open(f"analysis/bubble/{graphtype}_{mutype}_sv.tsv")
    lines = infile.readlines()

    merger = PdfFileMerger(strict=False)
    if len(lines) > 501:
        selsamp = sample(range(0, len(lines)), 500)
    else:
        selsamp = range(0, len(lines))

    for ind, line in enumerate(lines):
        if ind in selsamp:
            line_comp = line.strip().split()
            # if biallelic
            # 1 165873 AltDel 497 2 s1 s2 s120535 s3
            if mutype == "biallelic":
                start_node = line_comp[5]
                stop_node = line_comp[-1]
            # 1_575817        9       440     AltIns  s46,s145319,s48
            else:
                start_node = line_comp[-1].split(",")[0]
                stop_node = line_comp[-1].split(",")[-1]

            graph_viz.visualize_graph(graphtype, graphcomp, graph, nodeinf, start=start_node, stop=stop_node)
            merger.append(f"{graphtype}_{start_node}_{stop_node}.pdf")
            os.remove(f"{graphtype}_{start_node}_{stop_node}")
            os.remove(f"{graphtype}_{start_node}_{stop_node}.pdf")

    merger.write(f"analysis/bubble/{graphtype}_{mutype}_sv_viz.pdf")


if __name__ == "__main__":
    main()
