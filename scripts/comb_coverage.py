#!/usr/bin/env python

"""
Combine remapping file into a matrix 


Output: Matrix with format as follow

*Node coverage file*
|Node| Assemb1| Assemb2| Assemb3|
|s1| 0| 1| 0.125| 2|

*Edge coverage file*
|Parent | Child | Assemb1 | Assemb1|
|s1|s2|0|1|0|

Number indicates alignment coverage for given node or edge



"""

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Combine remapping file")
    parser.add_argument("-g", "--graph", help="Graph type")
    parser.add_argument("-a", "--assembly", help="Graph type", nargs="+")
    return parser.parse_args()


def parse_node_coverage(line):
    # S       s34     CGTGACT LN:i:7  SN:Z:1  SO:i:122101     SR:i:0  dc:f:0
    # "nodeid","nodelen","chromo","pos","rrank",assemb
    """
    Parse the gaf alignment
    Input: line from gaf alignment
    Output: tuple of nodeid, nodelen, start_chromo, start_pos, coverage

    """
    line_comp = line.strip().split()
    nodeid = line_comp[1]
    nodelen = len(line_comp[2])
    start_chromo = line_comp[4].split(":")[2]
    start_pos = line_comp[5].split(":")[2]
    rrank = line_comp[-2].split(":")[2]
    coverage = line_comp[-1].split(":")[2]

    return nodeid, nodelen, start_chromo, start_pos, rrank, coverage


def parse_edge_coverage(line):
    # L       s1      +       s133016 +       0M      SR:i:2  L1:i:165873     L2:i:2  dc:f:0
    linecomp = line.strip().split()
    if linecomp[2] == "-" or linecomp[4] == "-":
        parent = linecomp[3]
        child = linecomp[1]
    else:
        parent = linecomp[1]
        child = linecomp[3]
    cov_raw = int(linecomp[-1].split(":")[2])
    cov_rep = 1 if cov_raw > 0 else 0
    return (parent, child, cov_rep)


if __name__ == "__main__":
    args = parse_args()
    graph = args.graph
    assembly = args.assembly
    for assemb in assembly:
        infile = open(f"remap/{graph}/{assemb}_{graph}.gaf")

        # output node coverage
        if assemb == assembly[0]:
            combcov = pd.DataFrame([parse_node_coverage(line) for line in infile if line.startswith("S")],
                                   columns=["nodeid", "nodelen", "start_chromo", "start_pos", "rrank", assemb])
            infile.close()
        else:
            addcov = pd.DataFrame([[parse_node_coverage(line)[0], parse_node_coverage(line)[-1]] for line in infile if line.startswith("S")],
                                  columns=["nodeid", assemb])
            combcov = pd.merge(combcov, addcov, on=["nodeid"], how="outer")
            infile.close()

        # output edge coverage
        # Should take this out as a function, now just repetition
        infile = open(f"remap/{graph}/{assemb}_{graph}.gaf")

        if assemb == assembly[0]:
            combedge = pd.DataFrame([parse_edge_coverage(line) for line in infile if line.startswith("L")],
                                    columns=["parent_node", "child_node", assemb])
            infile.close()
        else:
            addedge = pd.DataFrame([parse_edge_coverage(line) for line in infile if line.startswith("L")],
                                   columns=["parent_node", "child_node", assemb])
            combedge = pd.merge(combedge, addedge, on=["parent_node", "child_node"], how="outer")
            infile.close()

    # write the node file
    combcov.fillna(0)
    combcov.to_csv(f"remap/{graph}_coverage.tsv",
                   sep=" ",
                   index=False)

    # write the edge file
    combedge.to_csv(f"remap/{graph}_edge_use.tsv",
                    sep=" ",
                    index=False)
