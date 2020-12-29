#!/usr/bin/env python

"""

Extract the biallelic SV from bubble in the graphs
and output the mutation types

"""
import argparse
from dataclasses import dataclass


@dataclass
class svid:
    """

    Object to store the sv information

    """

    chromo: str
    pos: int
    svtype: str
    reflen: int
    nonreflen: int
    sourcenode: str
    bubref: str
    bubnonref: str
    sinknode: str


def sv_det(line):
    """
    Determine SV types 
    Input: Each line of gfatools bubble output
    Output: Insertions or Deletions SV based on comparison
            between ref and non-ref allele

    """
    # 1 348029 4 2 s1,s2,s135091,s3
    chromo, pos, nodes, _, nodelist = line.strip().split()
    bubble = nodelist.split(",")[1:-1]
    sourcenode = nodelist.split(",")[0]
    sinknode = nodelist.split(",")[-1]
    if nodes in ["3", "4"]:
        if nodes == "3":
            rrank, nodelen = nodeinf[bubble[0]]
            if rrank > 0:
                svtype = "Insertion"
                reflen = 0
                nonreflen = nodelen
                bubnonref = bubble[0]
                bubref = 0
            else:
                svtype = "Deletion"
                reflen = nodelen
                nonreflen = 0
                bubnonref = 0
                bubref = bubble[0]
        elif nodes == "4":
            # do not consider looping in the bubble
            if bubble[0] == bubble[1]:
                return None
            for bub in bubble:
                rrank, nodelen = nodeinf[bub]
                if rrank > 0:
                    nonreflen = nodelen
                    bubnonref = bub
                else:
                    reflen = nodelen
                    bubref = bub
            svtype = "AltDel" if nonreflen < reflen else "AltIns"
        return chromo, pos, svtype, reflen, nonreflen, sourcenode, bubref, bubnonref, sinknode
    else:
        return None


if __name__ == "__main__":

    # parse the assembly
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--assembly", help="Assembly to process")
    args = parser.parse_args()
    assembly = args.assembly

    # container of the node including len and rank
    nodeinf = {}

    # s2 1389 1 348029 0
    with open(f"graph/{assembly}_graph_len.tsv") as infile:
        for line in infile:
            nodeid, nodelen, chromo, pos, rrank = line.strip().split()
            nodeinf[nodeid] = [int(rrank), int(nodelen)]

    with open(f"analysis/bubble/{assembly}_biallelic_bubble.tsv") as infile:
        for line in infile:
            if sv_det(line):
                print(*sv_det(line))
