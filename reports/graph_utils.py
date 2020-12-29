#!/usr/bin/env python
from collections import defaultdict


def revstrand(strand):
    """
    return + if strand - and vice versa

    """
    if strand == "+":
        return "_"
    elif strand == "-":
        return "+"


def revcomp(seq):
    basedict = {"A": "T", "G": "C", "T": "A", "C": "G"}
    seq = seq.upper()
    return "".join(basedict.get(x, x) for x in reversed(seq))


def parse_graph(input):
    grseq = {}
    graph = defaultdict(dict)

    # with open(f"graph/{assembly}_graph.gfa") as infile:
    with open(input) as infile:
        for line in infile:
            # parse the node sequences
            if line.startswith("S"):
                line_comp = line.strip().split()
                grseq[line_comp[1]] = line_comp[2]
            # parse the edge information
            # L       s1      +       s2      +       0M      SR:i:0  L1:i:426636     L2:i:961
            elif line.startswith("L"):
                parent, strand1, child, strand2 = line.strip().split()[1:5]
                if strand1 == "+" and strand2 == "+":
                    graph[parent][child] = [strand1, strand2]
                else:
                    graph[child][parent] = [
                        revstrand(strand2), revstrand(strand1)]
    return grseq, graph
