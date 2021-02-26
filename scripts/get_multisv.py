#!/usr/bin/env python
"""

Given minigraph multiallelic bubble 
Output the SV from the graphs including the mutation types (stdout)

"""

import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", "--assemb", help="assembly")
    return parser.parse_args()


def DFS(graph, start, end, path=None):
    if path is None:
        path = [start]
    else:
        path = path + [start]
    if start == end:
        paths.append(path)
    else:
        for node in graph[start]:
            if node not in path:
                DFS(graph, node, end, path)


@dataclass
class node():

    """

    class to store length, rank and colour from each node

    """
    rrank: int = 0
    nodelen: int = 0
    nodecol: str = ""


@dataclass
class bubbles():
    """

    class to store the path in the bubbles
    """

    # List length for each path
    pathlen: List[int] = field(default_factory=list)
    # List of rank List in path
    allrank: List[List] = field(default_factory=list)
    # bool whether path is ref path or not
    refpath: List[bool] = field(default_factory=list)
    # ref path length
    lenref: int = 0
    # List of rank List of colour in each path
    colour: List[List] = field(default_factory=list)
    # List of SV types
    svtype: List[str] = field(default_factory=list)


def get_sv(pathlen, allrank, refpath, lenref, colour):
    if refpath:
        return "ref"
    else:
        # Check if the path realistic
        # at least have one colour
        if colour:
            # if all nodes ref but not complete
            # it is deletions
            if all(x == 0 for x in allrank):
                return "Deletions"
            else:
                # Insertions when svlen longer than reflen
                if pathlen >= lenref:
                    if lenref == 0:
                        return "Insertions"
                    else:
                        # Altins when sv contains ref sequence
                        return "AltIns"
                else:
                    return "AltDel"

        else:
            return "not_path"


if __name__ == "__main__":
    args = parse_args()
    assemb = args.assemb

    # dict to get node list, length, rank and colour
    nodeinf = {}

    # s2 1389 1 348029 0
    with open(f"graph/{assemb}_graph_len.tsv") as infile:
        for line in infile:
            nodeid, nodelen, chromo, pos, rrank = line.strip().split()
            # nodeinf[nodeid] = [int(rrank), int(nodelen)]
            nodeinf[nodeid] = node(rrank=int(rrank), nodelen=int(nodelen))

    with open(f"analysis/colour_node/{assemb}_nodecol.tsv") as infile:
        # skip header
        next(infile)
        for line in infile:
            nodeid, nodecol = line.strip().split()
            nodeinf[nodeid].nodecol = nodecol

    # parse edges in the graph
    graph = defaultdict(list)
    with open(f"graph/{assemb}_graph.gfa") as infile:
        for line in infile:
            if line.startswith("L"):
                parent, strand1, child, strand2 = line.strip().split()[1:5]
                if strand1 == "+" and strand2 == "+":
                    graph[parent].append(child)
                else:
                    graph[child].append(parent)

    # determine types for each multiallelic bubble
    with open(f"analysis/bubble/{assemb}_multiallelic_bubble.tsv") as infile:
        for line in infile:
            chromo, pos, *_ = line.strip().split()
            nodes = line.strip().split()[4].split(",")
            start, *_, stop = nodes
            paths = []
            # do DFS from start to stop and
            # store all enumerated paths in paths
            DFS(graph, start, stop)
            source_node = int(paths[0][0].replace("s", ""))
            sink_node = int(paths[0][-1].replace("s", ""))
            # container to store properties from each path in paths
            bubble = bubbles()
            # Checking propertis from each path in paths
            for path in paths:
                totlen = 0
                rank_list = []
                colour_list = []
                for ind, comp in enumerate(path):
                    # determine the path length
                    # pathlen exclude source and sink
                    if ind > 0 and ind < len(path) - 1:
                        totlen = totlen + nodeinf[comp].nodelen
                    # determine which colour consistent across nodes
                    colour_list.append(set(nodeinf[comp].nodecol.split(",")))
                    # store the rrank list
                    rank_list.append(nodeinf[comp].rrank)
                bubble.colour.append(set.intersection(*colour_list))
                bubble.pathlen.append(totlen)
                bubble.allrank.append(rank_list)
                # determine which is the ref path and its length
                if all(x == 0 for x in rank_list) and len(path) == (sink_node - source_node + 1):
                    bubble.refpath.append(True)
                    bubble.lenref = totlen
                else:
                    bubble.refpath.append(False)

            # Determine sv type
            for ind, path in enumerate(paths):
                bubble.svtype.append(
                    get_sv(pathlen=bubble.pathlen[ind],
                           allrank=bubble.allrank[ind],
                           refpath=bubble.refpath[ind],
                           lenref=bubble.lenref,
                           colour=bubble.colour[ind])
                )

            for i in range(len(paths)):
                # not considering not_path as sv
                if bubble.svtype[i] != "not_path" and bubble.svtype[i] != "ref":
                    outvar = (
                        f"{chromo}_{pos}\t{bubble.lenref}\t"
                        f"{bubble.pathlen[i]}\t{bubble.svtype[i]}\t"
                        f"{','.join(paths[i])}"
                    )

                    print(outvar)
