#!/usr/bin/env python
"""
Intersect sv breakpoint with GFF annotations 
and report the most specific features from each breakpoints

intergenic --> intronic ---> exonic ---> CDS

"""
import pybedtools
import re
from collections import defaultdict


def get_feature_table(mergedbreak):
    """
    parse the gff line to get name and ID from each feature
    """
    for feature in mergedbreak:
        annot = [x.split("=") for x in feature.fields[-1].split(";")]
        try:
            annot2 = {x1: x2 for x1, x2 in annot}
            ID = annot2.get("ID", "noid")
            if ID not in "noid":
                ID = ID.split(":")[1]
            Name = annot2.get("Name", "noname")
            yield [feature[3], feature[6], ID, Name]
        except:
            continue


def extract_important_feature(prevfeat, curfeat):
    """

    Return the most important feature from gff
    """
    priority = {"CDS": 4, "exon": 3, "mRNA": 2, "gene": 1}
    prevprior = priority.get(prevfeat, 0)
    curprior = priority.get(curfeat, 0)

    if curprior == 0 and prevprior == 0:
        return "intergenic"
    if curprior >= prevprior:
        return curfeat
    else:
        return prevfeat


def get_feat_sel(svtype, featid, featname, featsel):
    """
        Extract gene name from gff features 
    """
    if re.search(r"gene", svtype):
        # report gene name which likely to be the shortest
        featsel = featid if len(featname) > len(featid) else featname
        # if no name, report the gene ID
        if featsel == "noname":
            featsel = featid
    return featsel


if __name__ == "__main__":
    graph = snakemake.wildcards["asb"]
    breakpointfile = pybedtools.BedTool(snakemake.input[0])
    annotfile = pybedtools.BedTool(snakemake.params["gffinput"])
    mergedbreak = breakpointfile.intersect(annotfile, wb=True, loj=True, stream=True)
    breakannot = {}
    combfeat = get_feature_table(mergedbreak)

    # process first line
    firstcomp = next(combfeat)
    svid, svtype, featid, featname = firstcomp
    typesel = extract_important_feature("", svtype)
    idprev = svid
    featsel = featid if len(featname) > len(featid) else featname
    featsel = get_feat_sel(svtype, featid, featname, featsel)

    # iterate on all features
    for comp in combfeat:
        svid, svtype, featid, featname = comp
        idnow = svid
        # get the most important features
        # if the svid the same
        if idnow == idprev:
            typesel = extract_important_feature(typesel, svtype)
            featsel = get_feat_sel(svtype, featid, featname, featsel)
        else:
            breakannot[idprev] = [typesel, featsel]
            # reset
            typesel = extract_important_feature("", svtype)
            featsel = featid if len(featname) > len(featid) else featname
            featsel = get_feat_sel(svtype, featid, featname, featsel)
        idprev = idnow
    # add the last component
    breakannot[idprev] = [typesel, featsel]
    outname = snakemake.output[0]
    with open(outname, "a") as outfile:
        for key, value in breakannot.items():
            print(key, *value, file=outfile)
