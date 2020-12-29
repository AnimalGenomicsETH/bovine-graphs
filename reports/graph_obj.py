#!/usr/bin/env python
"""

Class of the directed graph produced by minigraph

Consist of three classess:

    Graph : collection of nodes and edges
    Node  : node with attributes
    Edge  : edge with attributes

"""

from dataclasses import dataclass
from collections import defaultdict
from graph_utils import revcomp


@ dataclass
class Node():
    """

    Store node properties

    Attributes:

    nodeid : node identifier (s1,s2..)
    nodelen : length of sequences in the node (bp)
    nodeseq : nucleotide sequences in the node
    nodechr : chromosomal location node
    nodecoord : coordinate of the node
    nodecol : node label
    noderank : order of inclusion of node (0 reference)
    """

    nodeid: str = ""
    nodelen: int = 0
    nodeseq: str = ""
    nodechr: str = 0
    nodecoord: int = 0
    nodecol: str = ""
    noderank: int = 0

    @ classmethod
    def from_rgfa(cls, line, include_seq=1):
        """
        Alternative constructor from rGFA output
        Input: node line from rgfa, include_seq: 0 will not include seq in the graph
        Output: Node object

        """
        # S       s15     A       LN:i:1  SN:Z:1  SO:i:416249     SR:i:0
        line_comp = line.strip().split()
        nodeid = line_comp[1]
        if include_seq:
            nodeseq = line_comp[2]
        else:
            nodeseq = ""
        nodelen, nodechr, nodecoord, noderank = [x.split(":")[-1] for x in line_comp[3:]]
        nodecol = ""
        return cls(nodeid=nodeid,
                   nodeseq=nodeseq,
                   nodelen=int(nodelen),
                   nodechr=nodechr,
                   nodecoord=int(nodecoord),
                   nodecol=nodecol,
                   noderank=int(noderank))


@ dataclass
class Edge():
    """

    Store edge properties
    Edges always directed

    parent ---> child

    Strand denotes orientation (+/-)

    """

    parent: Node
    child: Node
    strand_parent: str
    strand_child: str

    @ classmethod
    def from_rgfa(cls, line, nodes):
        # L       s1      +       s2      +
        link, node1, strand1, node2, strand2, *_ = line.strip().split()
        node1 = int(node1[1:])
        node2 = int(node2[1:])
        # Only work on complete rGFA
        # assuming that nodes added with order preserved
        if "-" in strand1 or "-" in strand2:
            # need to reverse
            parent = nodes[node2 - 1]
            child = nodes[node1 - 1]
            strand1 = revcomp(strand2)
            strand2 = revcomp(strand1)
        else:
            parent = nodes[node1 - 1]
            child = nodes[node2 - 1]
        return cls(parent=parent, child=child, strand_parent=strand1, strand_child=strand2)


class Graph():
    """

    Store graph properties

    Attributes
    nodes : collection of nodes object
    edges: collection of edges object

    Methods
    add_node : add new node in the graph
    add_edge: add new edge in the graph
    graph["s12"]: slice the graph at given node
    graph_len: total pangenome length
    nonref_len: total non-ref sequences in the graph
    conv_adjlist: convert edges into an adjacency list representation



    """

    def __init__(self, nodes=None, edges=None):
        if nodes:
            self.nodes = nodes
        else:
            self.nodes = []

        if edges:
            self.edges = edges
        else:
            self.edges = []

    def add_nodes(self, node):
        """

        add a node to the graph

        Params:
        A node object

        Returns: None, but will add node to the nodes

        """

        self.nodes.append(node)

    def add_edges(self, edge):
        """

        add an edge to graph

        Params: an edges object

        Returns: None, but will add edge to the edges
        """

        self.edges.append(edge)

    def __getitem__(self, nodeid):
        """
        Slice node object for particular ID e.g.  Graph["s12"]
        """

        # convert id to int and substract 1 to reflect the ordering
        node_order = int(nodeid[1:]) - 1
        return self.nodes[node_order]

    def __str__(self):
        return (f"Graph contains {len(self.nodes)} nodes and {len(self.edges)} edges \n"
                f"Total graph length is {self.graph_len} bp with {self.nonref_len} non-reference bases")

    @ property
    def graph_len(self):
        """

        Calculate the length of graph

        Params: None
        Output: length of graphs (bp)

        """
        totlen = 0
        for node in self.nodes:
            totlen += node.nodelen
        return totlen

    @ property
    def nonref_len(self):
        """
        Calculate length of the non-ref sequence in the graph

        Params: None
        Output: Length non-ref sequences (bp)

        """
        nonref = 0
        for node in self.nodes:
            if node.noderank > 0:
                nonref += node.nodelen
        return nonref

    def conv_adjlist(self):
        """

        Convert edges representation into a adjacency list

        Params: None
        Output: Dict with parent as key and list containing all childs

        """

        adj_list = defaultdict(list)

        for edge in self.edges:
            adj_list[edge.parent].append(edge.child)

    @classmethod
    def construct_from_rgfa(cls, graphfile, Node=Node, Edge=Edge, include_seq=0):
        """ construct graph object from Minigraph output.

        :Input: graph: rgfa output from Minigraph
        :Input: Node : node class
        :Input: Edges : edge class
        :Input: include_seq : 1 if want to store sequences in nodes, default 0 (not storing sequences)
        :returns: Graph object

        """

        with open(graphfile) as infile:
            nodes_list = []
            edges_list = []
            for line in infile:
                # add nodes to the graph
                if line.startswith("S"):
                    node = Node.from_rgfa(line, include_seq)
                    nodes_list.append(node)
                # add edges to the graph
                if line.startswith("L"):
                    edge = Edge.from_rgfa(line=line, nodes=nodes_list)
                    edges_list.append(edge)
        return cls(nodes=nodes_list, edges=edges_list)
