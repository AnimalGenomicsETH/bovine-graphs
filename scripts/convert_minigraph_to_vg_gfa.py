#!/usr/bin/env python

from graph_obj import Node, Edge, Graph
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", help="input file")
    parser.add_argument("-o", "--output", help="output file")
    return parser.parse_args()


def create_gfa_node(graph):
    """ Create gfa input compatible with vg

    :input: graph object
    :return: node with numeric nodeid,  updated path

    """

    global pathlist

    for node in graph.nodes:
        if node.noderank == 0:
            chrnow = node.nodechr
            nodeid = node.nodeid.replace("s", "")
            pathlist[chrnow] = ",".join([pathlist[chrnow], f"{nodeid}+"])
        yield(f"S\t{nodeid}\t{node.nodeseq}\n")


def create_gfa_edges(graph):
    """

    Create a simplified gfa input for vg
    L parent strand child strand

    """

    for edge in graph.edges:
        parent_num = edge.parent.nodeid.replace("s", "")
        child_num = edge.child.nodeid.replace("s", "")
        yield((f"L\t{parent_num}\t{edge.strand_parent}\t"
               f"{child_num}\t{edge.strand_child}\n"))


if __name__ == "__main__":
    args = parse_args()
    input_graph = args.input
    output_graph = args.output
    graph = Graph.construct_from_rgfa(input_graph, include_seq=1)

    pathlist = defaultdict(str)

    with open(output_graph, "a") as outfile:
        # wrote node and updated pathlist
        for comp in create_gfa_node(graph):
            outfile.write(comp)
        # wrote edges
        for comp in create_gfa_edges(graph):
            outfile.write(comp)
        # wrote path
        for pathid, path in pathlist.items():
            print(f"P\t{pathid}\t{path[1:]}\t*", file=outfile)
