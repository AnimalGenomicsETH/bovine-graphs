#!/usr/bin/env python
"""

Scripts to report the graph statistics

"""
from graph_obj import Graph, Edge, Node
from dataclasses import dataclass
from collections import defaultdict
from graph_utils import revcomp


def node_classifier(graph):
    """
    Given the graph it will output the node augmented from other assembly
    Input: Graph object
    Output: Dict with rank as key and the tuple of node count and seq len as items

    """
    node_counter = defaultdict(lambda: [0, 0])
    for node in graph.nodes:
        noderank = node.noderank
        node_counter[noderank][0] += 1
        node_counter[noderank][1] += node.nodelen
    return node_counter


def edge_classifier(graph):
    """
    Given graph object, it will output the edges type
    Input: Graph object
    Output: List the count of ref-ref, ref-non-ref, and non-ref-non-ref edges

    """
    edge_counter = {"RefRef": 0, "RefNonref": 0, "NonrefNonref": 0}
    for edge in graph.edges:
        all_ranks = [edge.parent.noderank, edge.child.noderank]
        if all([x == 0 for x in all_ranks]):
            edge_counter["RefRef"] += 1
        elif all([x > 0 for x in all_ranks]):
            edge_counter["NonrefNonref"] += 1
        else:
            edge_counter["RefNonref"] += 1
    return edge_counter


def graph_stat_report(assembly, assemb_comp):
    graphfile = f"graph/{assembly}_graph.gfa"
    combgraph = Graph.construct_from_rgfa(graphfile=graphfile)
    report_string = "\n\n\n### *Graph statistics*\n\n\n"
    report_string += ("| Graph parameters | Count | Length (bp) | \n"
                      "|----|----|----| \n")
    # Add overall graph len
    report_string += f"| All nodes | {len(combgraph.nodes):,} | {combgraph.graph_len:,} | \n"
    all_nodes = node_classifier(combgraph)
    # Add ref len
    ref_count = all_nodes[0][0]
    ref_len = all_nodes[0][1]
    report_string += f"| Reference nodes| {ref_count:,} | {ref_len:,}| \n "
    nonref_count = sum([values[0] for key, values in all_nodes.items() if key > 0])
    nonref_len = sum([values[1] for key, values in all_nodes.items() if key > 0])
    # Add non-reference graph len
    report_string += f"| Non-reference nodes| {nonref_count:,} | {nonref_len:,}| \n "
    # Add from non-ref from each assembly
    for ind, [key, values] in enumerate(all_nodes.items()):
        if key != 0:
            report_string += f"|Added from {assemb_comp[ind]}|{int(values[0]):,} | {int(values[1]):,}| \n"
    # Add total number of edges
    total_edges = sum(edge_classifier(combgraph).values())
    report_string += f"|Total Edges| {total_edges} | ratio: {total_edges/len(combgraph.nodes):.4f}|\n"
    # Add edge types
    for key, value in edge_classifier(combgraph).items():
        report_string += f"|Edge {key}|{value:,}|0| \n"
    return report_string
