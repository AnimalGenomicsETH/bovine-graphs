#!/usr/bin/env python

from graph_obj import Graph, Node, Edge
from graph_utils import revcomp
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--bialfile", help="biallelic sv stat file")
    parser.add_argument("-m", "--multifile", help="multiallelic sv stat file")
    parser.add_argument("-g", "--graphfile", help="graphfile")
    parser.add_argument("-o", "--outputfile", help="output file")
    return parser.parse_args()


def parse_bial_sv(stat, graph):
    nodeid = stat.strip().split()[3].split(",")[1]
    strand = stat.strip().split()[4].split(",")[1]
    sv_name = stat.strip().split()[0]
    if strand == "-":
        sequence = revcomp(graph[nodeid].nodeseq)
    else:
        sequence = graph[nodeid].nodeseq
    return (sv_name, sequence)


def parse_multi_sv(stat, graph):
    linecomp = stat.strip().split()
    sv_name = linecomp[0]
    nodelist = linecomp[3].split(",")[1:-1]
    strandlist = linecomp[4].split(",")[1:-1]
    allcomp = zip(nodelist, strandlist)
    sequence = ""
    for nodeid, strand in allcomp:
        if strand == "-":
            sequence += revcomp(graph[nodeid].nodeseq)
        else:
            sequence += graph[nodeid].nodeseq
    return (sv_name, sequence)


if __name__ == "__main__":
    args = parse_args()
    bialfile = args.bialfile
    multifile = args.multifile
    graphfile = args.graphfile
    outputfile = args.outputfile
    # construct graph
    graph = Graph.construct_from_rgfa(graphfile, include_seq=1)

    # process biallelic sv
    with open(bialfile) as infile, open(outputfile, "a") as outfile:
        for line in infile:
            sv_name, sequence = parse_bial_sv(line, graph)
            outfile.write(f">{sv_name}\n{sequence}\n")
    # process multiallelic sv
    with open(multifile) as infile, open(outputfile, "a") as outfile:
        for line in infile:
            sv_name, sequence = parse_multi_sv(line, graph)
            outfile.write(f">{sv_name}\n{sequence}\n")
