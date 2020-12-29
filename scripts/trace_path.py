#!/usr/bin/env python
"""

For paths in the bubles,
output which assemblies touching path

Input: Edge coverage file
Output: stdout each paths and assembly it derived

"""

from collections import defaultdict
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-g", "--graph", help="graph to calculate")
    parser.add_argument("-a", "--assembly", help="assemblies in the graphs", nargs="+")
    return parser.parse_args()

def parse_edge_coverage(covfile):
    """

    Create nested dict with [parent][child][assembly] from coverage file

    """
    basedict = lambda: defaultdict(basedict)
    edge_cover = basedict()

    with open(covfile) as infile:
        _, _, *assemb = next(infile).strip().split()

        for line in infile:
            parent, child, *covers = line.strip().split()

            for ind, cover in enumerate(covers):
                edge_cover[parent][child][assemb[ind]] = int(cover)

    return edge_cover


def generate_path(start, ref, nonref, stop):
    """

    Given source, sink, and ref/non-ref nodes enumerate all possible paths

    """
    ref_path = [x for x in [start, ref, stop] if x != "0"]
    nonref_path = [x for x in [start, nonref, stop] if x != "0"]
    return [ref_path, nonref_path]


def path_traverse(path, anim, edge_cover):
    """

    Check whether path is traversed for a given assembly

    """

    all_cover = []
    for ind, node in enumerate(path[:-1]):
        parent = node
        child = path[ind + 1]
        cover = edge_cover[parent][child][anim]
        all_cover.append(cover)

    if all(all_cover):
        return anim
    else:
        return ""


def main():
    anims = ["UCD", "Angus", "Highland", "OBV", "Brahman", "Yak"]

    args = parse_args()
    anims = args.assembly
    graph = args.graph

   # dict of the edge coverage
    edge_cover = parse_edge_coverage(f"remap/{graph}_edge_use.tsv")

    def traverse_all_anims(ref_path, nonref_path, ref_len, nonref_len, mutype, anim=anims, svmode="biallelic"):
        """
        Going to the each assembly and check whether it traverse the ref and non-ref path

        """
        # iterate for each assembly whether its support ref or nonref
        ref_list = [path_traverse(ref_path, anim, edge_cover=edge_cover) for anim in anims]
        nonref_list = [path_traverse(nonref_path, anim, edge_cover=edge_cover) for anim in anims]

        # if no assembly can traverse just write noassemb
        if all(x == "" for x in nonref_list):
            nonref_list = ["noassemb"]

        if all(x == "" for x in ref_list):
            ref_list = ["noassemb"]

        # return with svid, ref paths and the label
        print(svmode, "\t", svid, "\t", mutype, "\t", ref_len, "\t", nonref_len, "\t",
              ",".join(ref_path), "\t", ",".join([x for x in ref_list if x]), "\t",
              ",".join(nonref_path), "\t", ",".join([x for x in nonref_list if x]))

    # process biallelic bubble
    with open(f"analysis/bubble/{graph}_biallelic_sv.tsv") as infile:
        for line in infile:
            # 1 165873 AltDel 497 2 s1 s2 s133016 s3
            linecomp = line.strip().split()
            svid = "_".join(linecomp[:2])
            mutype, ref_len, nonref_len = linecomp[2:5]
            start, ref, nonref, stop = linecomp[5:]
            ref_path, nonref_path = generate_path(start, ref, nonref, stop)

            traverse_all_anims(ref_path, nonref_path, ref_len, nonref_len, mutype)

    # process multiallelic bubble
    with open(f"analysis/bubble/{graph}_multiallelic_sv.tsv") as infile:
        for line in infile:
            # 1_535561        1898    5138    AltIns  s35,s123493,s38
            svid, ref_len, nonref_len, mutype, nonref_path = line.strip().split()
            # create ref path
            nonref_path = nonref_path.split(",")
            start_node, *_, stop_node = nonref_path
            startid = int(start_node[1:])
            stopid = int(stop_node[1:])
            ref_path = [start_node, *[f"s{x}"for x in range(startid + 1, stopid)], stop_node]

            traverse_all_anims(ref_path, nonref_path, ref_len, nonref_len, mutype, svmode="multiallelic")


if __name__ == "__main__":
    main()
