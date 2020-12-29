#!/usr/bin/env python
""""

Given biallelic bubble collections
output the sequences with correct sequences orientation
e.g if (-) reverse complemented

output:

100 bp (upstream, source) --- bubble sequences --- 100 bp (downstream, sink)

"""

import argparse
from collections import defaultdict
from graph_utils import revstrand, revcomp, parse_graph


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--assembly", help="Assembly to process")
    return parser.parse_args()


def get_svseq(source, norefbub, sink):
    """

    Given the sequence from source and bubble and sink
    stitch them together
    Return: the concatenated sequences

    """
    chainstrand = graph[source][norefbub] + [graph[norefbub][sink][0]]
    # use the first 100 bp, if the ref node less than use all refnode
    try:
        totseq = grseq[source][-100:] if chainstrand[0] == "+" else revcomp(
            grseq[source])[-100:]
    except:
        totseq = grseq[source] if chainstrand[0] == "+" else revcomp(
            grseq[source])

    totseq += grseq[norefbub] if chainstrand[1] == "+" else revcomp(
        grseq[norefbub])

    try:
        totseq += grseq[sink][:100] if chainstrand[2] == "+" else revcomp(grseq[sink])[
            :100]

    except:
        totseq += grseq[sink] if chainstrand[2] == "+" else revcomp(
            grseq[sink])

    return totseq


if __name__ == "__main__":
    args = parse_args()
    assembly = args.assembly
    # container to store node and edge sequences
    grseq, graph = parse_graph(f"graph/{assembly}_graph.gfa")
    # with open(f"graph/{assembly}_graph.gfa") as infile:
    # for line in infile:
    # # parse the node sequences
    # if line.startswith("S"):
    # line_comp = line.strip().split()
    # grseq[line_comp[1]] = line_comp[2]
    # # parse the edge information
    # #L       s1      +       s2      +       0M      SR:i:0  L1:i:426636     L2:i:961
    # elif line.startswith("L"):
    # parent, strand1, child, strand2 = line.strip().split()[1:5]
    # if strand1 == "+" and strand2 == "+":
    # graph[parent][child] = [strand1, strand2]
    # else:
    # graph[child][parent] = [revstrand(strand2), revstrand(strand1)]

    regfile = open(f"analysis/bubble/{assembly}_bialsv_stat.tsv", "a")
    errfile = open(f"analysis/bubble/{assembly}_bialsv_conflict.tsv", "a")
    bialseq = open(f"analysis/bubble/{assembly}_bialsv_seq.fa", "a")
    conseq = 0

    with open(f"analysis/bubble/{assembly}_biallelic_sv.tsv") as infile:
        for line in infile:
            line_comp = line.strip().split()
            chromo, pos, svtype, reflen, nonreflen, source, bubref, norefbub, sink = line_comp
            # 1 348029 AltDel 1389 2 s2 s1 s135091 s3
            # only consider non-deletion with length more than 100 bp
            if svtype != "Deletion" and int(nonreflen) >= 100:
                # check that there's no left and right link conflict
                if graph[source][norefbub][-1] == graph[norefbub][sink][0]:
                    totseq = get_svseq(source, norefbub, sink)
                    chainstrand = graph[source][norefbub] + \
                        [graph[norefbub][sink][0]]
                    conseq += 1
                    sv_name = f"b{conseq}_{chromo}_{pos}"
                    bialseq.write(f">{sv_name}\n")
                    bialseq.write(f"{totseq}\n")
                    regfile.write((f"{sv_name}\t{svtype}\t{len(totseq)}\t"
                                   f"{','.join([source,norefbub,sink])}\t"
                                   f"{','.join(chainstrand)}\n"))
                else:
                    print(chromo, pos, svtype,
                          reflen, nonreflen,
                          source, bubref, norefbub, sink,
                          graph[source][norefbub], graph[norefbub][sink],
                          file=errfile)

    regfile.close()
    errfile.close()
    bialseq.close()
