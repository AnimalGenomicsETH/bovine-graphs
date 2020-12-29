#!/usr/bin/env python


from graph_obj import Edge, Node, Graph
from collections import defaultdict
import argparse
import re


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--coverage", help="list of coverage file", nargs="+")
    parser.add_argument("-g", "--graph", help="graphfile")
    parser.add_argument("-o", "--output", help="output file")
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


def combine_node_coverage(nodefile):
    """

    Update nodecov dict with coverage file
    """
    global nodecov

    with open(nodefile) as infile:
        for line in infile:
            nodeid, coverage = line.strip().split()
            nodecov[nodeid] += int(coverage)


def combine_edge_coverage(edgefile):
    """

    Update edge file with coverage file
    """

    global edgecov

    with open(edgefile) as infile:
        for line in infile:
            parent, child, coverage = line.strip().split()
            edgecov[parent][child] += int(coverage)


if __name__ == "__main__":
    args = parse_args()
    coverage_file = args.coverage
    graphfile = args.graph
    output = args.output

    # dict of the node and edge
    nodecov, edgecov = parse_graph(graphfile)

    # combine node or edge files
    if re.search(r"node", coverage_file[0]):
        for covfile in coverage_file:
            combine_node_coverage(covfile)

        # write combined node coverage
        with open(output, "a") as outfile:
            for key, value in nodecov.items():
                outfile.write(f"{key}\t{value}\n")

    else:
        for covfile in coverage_file:
            combine_edge_coverage(covfile)

        # write combined edges coverage
        with open(output, "a") as outfile:
            for parent, child in edgecov.items():
                for key, value in child.items():
                    outfile.write(f"{parent} {key} {value}\n")
