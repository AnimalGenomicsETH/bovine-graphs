#!/usr/bin/env python

from graphviz import Digraph
from collections import defaultdict
from random import sample
import itertools
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--graph", help="graph to process")
    parser.add_argument("-c", "--graphcomp", help="assembly in the graph", nargs="+")
    parser.add_argument("--start", help="start/source node", type=str)
    parser.add_argument("--stop", help="stop/sink node", type=str)
    return parser.parse_args()


def generate_edges(ingraph):
    # databose on graph edges
    graph = defaultdict(list)
    with open(f"graph/{ingraph}_graph_link.tsv") as infile:
        for line in infile:
            if line.startswith("L"):
                link, parent, strand1, child, strand2, *_ = line.strip().split()
                if strand1 == "+" and strand2 == "+":
                    graph[parent].append(child)
                else:
                    graph[child].append(parent)
    return graph


def graph_info(graphlen):
    # database on nodes and its length
    nodeinf = {}

    # s2 1389 1 348029 0
    with open(f"graph/{graphlen}_graph_len.tsv") as infile:
        for line in infile:
            nodeid, nodelen, chromo, pos, rrank = line.strip().split()
            nodeinf[nodeid] = [int(rrank), int(nodelen), chromo, int(pos)]
    return nodeinf


def generate_colour(graphcomp):
    colist = ["lightsalmon1", "orange",
              "plum2", "palegreen2",
              "cyan2", "gold3",
              "firebrick2", "deepskyblue",
              "cadetblue3", "burlywood3"]
    colsel = []
    for ind, x in enumerate(itertools.cycle(colist)):
        if ind < len(graphcomp):
            colsel.append(x)
        else:
            break
    return colsel


def DFS(graph, start, end, path=None, paths=None, counter=0):
    """

    Depth first traversal from start to end node

    """

    if path is None:
        path = [start]
        paths = []
        counter += 1
    else:
        path = path + [start]
        counter += 1

    if start == end:
        paths.append(path)
    else:
        for node in graph[start]:
            if counter > 100:
                raise RuntimeError("Either start/stop not correct or recursion too deep")
            if node not in path:
                DFS(graph, node, end, path, paths, counter)
    return paths


def visualize_graph(graphtype, graphcomp, graph, nodeinf, start, stop, annot, outf="pdf"):
    """

    Visualize surrounding SVs from start-stop
    Input:
    Graphtype = name of graph, prefix of gfa
    Graphcomp = list of assembly in graph, List
    graph = graph object from generate_edges output
    nodeinf = node object from graph_info output
    start = start node to visualize "s1" 
    stop = stop node to visualize "s2"
    annot = annotation to add in the graph

    Output:
    None, but gives the figure of the SV
    """
    colsel = generate_colour(graphcomp)

    paths = DFS(graph, start, stop)

    # parse left and right

    # enumerate all nodes and edges from bubble

    all_node = list(itertools.chain(*paths))
    ref_node = list({x for x in all_node if nodeinf[x][0] == 0})
    alt_node = list({x for x in all_node if nodeinf[x][0] != 0})

    alt_node.sort(key=lambda x: int(x.replace("s", "")))
    ref_node.sort(key=lambda x: int(x.replace("s", "")))

    f = Digraph('bubble_multi', filename=f"{graphtype}_{start}_{stop}.dot", format=outf,
                graph_attr={'label': annot})
    f.attr(rankdir='LR', size='8,5')

    f.attr('node', shape="box")
    f.attr('node', color='bisque1', style="filled")
    for node in ref_node:
        f.node(node, label=f"{node}\n{nodeinf[node][1]} bp")

    for node in alt_node:
        f.attr('node', color=colsel[nodeinf[node][0]], style="filled")
        # f.node(node, label=f"{node}\n{nodeinf[node][1]} bp")
        f.node(node, label=f"{node}\n{nodeinf[node][1]} bp\n{nodeinf[node][2]}:{nodeinf[node][3]}")

    for parent in set(all_node):
        if parent not in [stop]:
            for child in graph[parent]:
                f.edge(parent, child)

    f.render(filename=f"{graphtype}_{start}_{stop}")


def main():
    args = parse_args()
    graphtype = args.graph
    graphcomp = args.graphcomp
    start = args.start
    stop = args.stop
    graph = generate_edges(graphtype)
    nodeinf = graph_info(graphtype)

    visualize_graph(graphtype, graphcomp, graph, nodeinf, start, stop, annot)


if __name__ == "__main__":
    main()
