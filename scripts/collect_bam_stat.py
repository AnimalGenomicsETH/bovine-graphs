#!/usr/bin/env python
"""
collect statistics from bam file

"""

import argparse
import pysam
import re


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--anims",
                        help="anims to call")
    parser.add_argument("-m", "--mode",
                        help="alignmentmode")
    return(parser.parse_args())


if __name__ == "__main__":
    args = parse_args()
    anims = args.anims
    mode = args.mode
    samfile = pysam.AlignmentFile(f"wgs/bam/{anims}_{mode}.bam", "rb")
    # dict to save the stat to count

    statcount = {"totmap": 0, "unmap": 0, "mapq0": 0, "mapq10": 0, "mapq60": 0,
                 "al99": 0, "alp": 0, "clip": 0, "suppl": 0, "secondary": 0,
                 "supl_tag": 0, "primary1": 0, "primary2": 0, "totun": 0}

    for read in samfile:
        statcount["totmap"] += 1
        if read.is_unmapped:
            statcount["unmap"] += 1
        else:
            if read.mapping_quality == 0:
                statcount["mapq0"] += 1
            if read.mapping_quality > 10:
                statcount["mapq10"] += 1
            if read.mapping_quality == 60:
                statcount["mapq60"] += 1
            edis = read.get_tag("NM")
            identity = (read.query_alignment_length - edis) / read.query_alignment_length
            if identity >= 0.99:
                statcount["al99"] += 1
            if re.search(r"S|H", read.cigarstring):
                statcount["clip"] += 1
            if not re.search(r"S|H", read.cigarstring) and identity == 1:
                statcount["alp"] += 1
            if read.is_supplementary:
                statcount["suppl"] += 1
            if read.is_secondary:
                statcount["secondary"] += 1
            try:
                read.get_tag("SA")
                statcount["supl_tag"] += 1
            except:
                try:
                    read.get_tag("XA")
                    if read.mapping_quality == 60:
                        statcount["primary2"] += 1
                except:
                    if read.mapping_quality > 10:
                        statcount["primary1"] += 1

    statcount["totun"] = statcount["primary1"] + statcount["primary2"]
    stat = ["totmap", "unmap", "mapq0", "mapq10", "mapq60", "al99",
            "alp", "clip", "suppl", "secondary", "supl_tag", "primary1", "primary2", "totun"]
    print(" ".join(str(statcount[x]) for x in stat), anims, mode)
