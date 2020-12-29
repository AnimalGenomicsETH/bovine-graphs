#!/usr/bin/env python
"""

Given multi-allelic bubble collection
output sequence by stiching sequences all nodes as a path
there should be no strand conflict

output:

100 bp upstream --- bubble path sequences --- 100 bp downstream (stdout)

"""

from collections import defaultdict
import argparse
from graph_utils import revstrand, revcomp, parse_graph


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--assembly", help="Assembly to process")
    return parser.parse_args()


def chain_node(graph, paths):
    """
    Given a multi node paths, 
    chain it together without strand conflict

    output: chained order e.g. +,-,- 
    where negative need to be reverse complemented

    """

    chainstrand = []
    for ind in range(1, len(paths) - 1):
        prevnode = paths[ind - 1]
        curnode = paths[ind]
        nextnode = paths[ind + 1]
        strandprev = graph[prevnode][curnode][-1]
        strandnext = graph[curnode][nextnode][0]
        if strandprev == strandnext:
            # create window prev, cur and next node
            # only stich when the strand consistents
            # if the paths only 3 nodes just append it directly (no sticthing)
            if len(paths) == 3:
                chainstrand.extend(graph[prevnode][curnode])
                chainstrand.extend(graph[curnode][nextnode][1])
            # when in the first window add the prev and cur strand
            elif ind == 1:
                chainstrand.extend(graph[prevnode][curnode])
            # when at the end add the current and next node
            elif ind == len(paths) - 2:
                chainstrand.extend(graph[curnode][nextnode])
            # otherwise add only the current node
            else:
                chainstrand.extend([graph[curnode][nextnode][0]])
        else:
            # report that there's conflict
            return chainstrand, strandprev, strandnext
        # report successful with the strand
    return chainstrand, strandprev, strandnext


def get_sv_seq(grseq, paths, chainstrand, ):
    """

    Given paths output the sequences 
    by concatenating with correct orientation

    Input:
    - grseq: dict map nodeid and sequences
    - paths: ordered nodes in a path
    - chainstrand: the corresponding list of orientation +/- 

    Return:
    Concatenated sequences


    """

    for ind, path in enumerate(paths):
        # add the 100 bp flanking ref for the start of the path
        if not ind:
            try:
                totseq = grseq[path][-100:] if chainstrand[0] == "+" else revcomp(
                    grseq[path])[-100:]
            except:
                totseq = grseq[path] if chainstrand[0] == "+" else revcomp(
                    grseq[path])
        else:
            # add 100 bp ref at the end of the path
            if ind == (len(paths) - 1):
                try:
                    totseq += grseq[path][:100] if chainstrand[2] == "+" else revcomp(grseq[path])[
                        :100]
                except:
                    totseq += grseq[path] if chainstrand[2] == "+" else revcomp(
                        grseq[path])
            else:
                totseq += grseq[path] if chainstrand[2] == "+" else revcomp(
                    grseq[path])
    return totseq


if __name__ == "__main__":
    args = parse_args()
    assembly = args.assembly
    grseq, graph = parse_graph(f"graph/{assembly}_graph.gfa")

    regfile = open(f"analysis/bubble/{assembly}_multisv_stat.tsv", "a")
    errfile = open(f"analysis/bubble/{assembly}_multisv_conflict.tsv", "a")
    multiseq = open(f"analysis/bubble/{assembly}_multisv_seq.fa", "a")

    with open(f"analysis/bubble/{assembly}_multiallelic_sv.tsv") as infile:
        sv_count = 0
        for line in infile:
            line_comp = line.strip().split()
            # 1_1579557       392     482     AltIns  s54,s55,s140149,s57
            pos, reflen, nonreflen, svtype, paths = line_comp
            paths = paths.split(",")
            if svtype != "Deletions" and int(nonreflen) > 100:
                # chain node with correct orientation
                chainstrand, strandprev, strandnext = chain_node(graph, paths)
                # output sequences if there's no conflict
                if len(chainstrand) == len(paths):
                    sv_count += 1
                    totseq = get_sv_seq(grseq, paths, chainstrand)
                    regfile.write(
                        f"m{sv_count}_{pos}\t{svtype}\t{len(totseq)}\t{','.join(paths)}\t{','.join(chainstrand)}\n")
                    multiseq.write(f">m{sv_count}_{pos}\n")
                    multiseq.write(f"{totseq}\n")
                # otherwise output conflict in errfile
                else:
                    print(*line_comp, chainstrand, strandprev,
                          strandnext, file=errfile)
