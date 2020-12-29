#!/usr/bin/env python

from graph_obj import Edge, Node, Graph
from collections import defaultdict
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--graph", help="graphfile")
    parser.add_argument("-a", "--alignment", help="alignment file")
    # output as list with first as node coverage file and second as the edge coverage file
    parser.add_argument("-o", "--outfile", help="output file", nargs="+")
    return parser.parse_args()


def parse_graph(graphfile):
    """

    Input GFA and output dict of nodes and edges

    """

    graph = Graph.construct_from_rgfa(graphfile)

    nodecov = {node.nodeid: 0 for node in graph.nodes}

    edgecov = defaultdict(dict)
    for edge in graph.edges:
        edgecov[edge.parent.nodeid][edge.child.nodeid] = 0

    return nodecov, edgecov


def coverage_updater(pathcomp):
    """

    Update node and edge coverage from reads

    """
    global nodecov
    global edgecov
    global added_node
    global added_edge

    for ind, comp in enumerate(pathcomp):
        if comp not in added_node:
            nodecov[comp] += 1
            added_node.append(comp)

            # if more than one nodes also calculate edge coverage
            # stop when at the end nodes
            if ind + 1 != len(pathcomp):
                parent = comp
                child = pathcomp[ind + 1]
                edge_id = f"{parent}_{child}"
                # only added edge never seen before in the same read

                if edge_id not in added_edge:
                    added_edge.append(edge_id)
                    try:
                        edgecov[parent][child] += 1
                    except:
                        try:
                            edgecov[child][parent] += 1
                        # only added edge which present in the initial graph
                        except:
                            continue


if __name__ == "__main__":
    args = parse_args()
    graphfile = args.graph
    infile = args.alignment
    outnode, outedges = args.outfile

    # dict of the node and edge
    nodecov, edgecov = parse_graph(graphfile)

    # only counting the unique reads
    # no doubling counting from the same reads
    readprev = ""
    added_node = []
    added_edge = []

    # process each line in gaf alignment
    with open(infile) as infile:
        for line in infile:
            read, qlen, qstart, qend, strand, pathlist, pathlen, pathstart, pathend, resid, alblock, mq, edis, dv, pid, cigar = line.strip().split()
            readnow = read
            pathcomp = pathlist.replace(">", "<").split("<")
            pathcomp = pathcomp[1:]

            if readnow == readprev:
                coverage_updater(pathcomp)

            else:
                # reset
                added_node = []
                added_edge = []
                # process the line

                coverage_updater(pathcomp)

        # reset read name
            readprev = readnow

    # output file
    with open(outnode, "a") as outfile:
        for key, value in nodecov.items():
            outfile.write(f"{key}\t{value}\n")

    with open(outedges, "a") as outfile:
        for parent, child in edgecov.items():
            for key, value in child.items():
                outfile.write(f"{parent} {key} {value}\n")
