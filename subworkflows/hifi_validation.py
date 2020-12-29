#!/usr/bin/env python
import glob
from pathlib import Path

parentdir = Path(srcdir("")).parents[0]
configfile: "config/config_hifi.yaml"
workdir: config["workdir"]

hifi_base = config["hifi_base"]
hifi_list = [x.split("/")[-1].split(".")[0] for x in glob.glob(hifi_base + "/*.fastq.gz")]
split_size = config["split_size"]
graph_list = config["graph_list"]

rule all:
    input:
        expand("validation/coverage/{hifi_anim}_{graph}_node_coverage.tsv", hifi_anim=hifi_list, graph=graph_list),
        expand("validation/coverage/{hifi_anim}_{graph}_edge_coverage.tsv", hifi_anim=hifi_list, graph=graph_list),
        expand("validation/coverage/{hifi_anim}_{graph}_bubble_support.tsv", hifi_anim=hifi_list, graph=graph_list)

checkpoint split_fasta:
    input:
        hifi_base + "/{hifi_anims}.fastq.gz"
    output:
        directory("validation/hifi_reads/split_{hifi_anims}")
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    params:
        split_size = split_size
    shell:
        """

        seqkit split2 -s {params.split_size} --threads {threads} -O {output} {input}


        """

rule map_hifi:
    input:
        fasta = "validation/hifi_reads/split_{hifi_anims}/{hifi_anims}.part_{part}.fastq.gz",
        graph = "graph/{graph}_graph.gfa"
    output:
        "validation/aligned/{hifi_anims}_{graph}/{hifi_anims}_{graph}_{part}.gaf"
    threads: 18
    resources:
        mem_mb = 3000,
        walltime = "04:00",
        disk_scratch = 20
    shell:
        """

        mkdir -p $TMPDIR/validation/hifi_reads/split_{wildcards.hifi_anims} && cp {input.fasta} $_
        mkdir -p $TMPDIR/graph && cp {input.graph} $_

        cd $TMPDIR

        graphaligner -x vg -t {threads} -g {input.graph} -f {input.fasta} -a out.gaf

        mv out.gaf $LS_SUBCWD/{output}


        """

rule calculate_coverage:
    input:
        alignment = "validation/aligned/{hifi_anims}_{graph}/{hifi_anims}_{graph}_{part}.gaf",
        graph = "graph/{graph}_graph.gfa"
    output:
        "validation/coverage/{hifi_anims}_{graph}/{hifi_anims}_{graph}_{part}_nodecov.tsv",
        "validation/coverage/{hifi_anims}_{graph}/{hifi_anims}_{graph}_{part}_edgecov.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    params:
        parentdir = parentdir
    shell:
        """

        {params.parentdir}/scripts/calculate_coverage_hifi.py -g {input.graph} -a {input.alignment} -o {output}

        """


def get_part(wildcards):
    """
    Return all splitted parts for a sample

    """
    selanim = wildcards.hifi_anims
    selgraph = wildcards.graph

    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    all_parts, = glob_wildcards(checkpoint_output + f"/{selanim}.part_{{part}}.fastq.gz")
    return [[f"validation/coverage/{selanim}_{selgraph}/{selanim}_{selgraph}_{part}_nodecov.tsv" for part in all_parts],
            [f"validation/coverage/{selanim}_{selgraph}/{selanim}_{selgraph}_{part}_edgecov.tsv" for part in all_parts]]


rule combined_node_hifi:
    input:
        node_list = lambda wildcards: get_part(wildcards)[0],
        graph = "graph/{graph}_graph.gfa"
    output:
        "validation/coverage/{hifi_anims}_{graph}_node_coverage.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "00:10"
    params:
        parentdir = parentdir
    shell:
        """

        {params.parentdir}/scripts/combine_hifi_coverage.py -g {input.graph} -c {input.node_list} -o {output}

        """

rule combined_edge_hifi:
    input:
        node_list = lambda wildcards: get_part(wildcards)[1],
        graph = "graph/{graph}_graph.gfa"
    output:
        "validation/coverage/{hifi_anims}_{graph}_edge_coverage.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "00:10"
    params:
        parentdir = parentdir
    shell:
        """

        {params.parentdir}/scripts/combine_hifi_coverage.py -g {input.graph} -c {input.node_list} -o {output}

        """

localrules: calculate_bubble_support
rule calculate_bubble_support:
    input:
        "analysis/bubble/{graph}_path_trace.tsv",
        "validation/coverage/{hifi_anims}_{graph}_edge_coverage.tsv"
    output: "validation/coverage/{hifi_anims}_{graph}_bubble_support.tsv"
    run:
        from collections import defaultdict

        def parse_edge_file(edgefile):
            combedge = defaultdict(dict)
            with open(edgefile) as infile:
                for line in infile:
                    parent, child, coverage = line.strip().split()
                    combedge[parent][child] = coverage
            return combedge

        def extract_coverage(combedge, paths):
            coverage = 0
            for ind, node in enumerate(paths[:-1]):
                parent = node
                child = paths[ind + 1]
                nodecover = int(combedge[parent][child])
                coverage = min([coverage, nodecover])
            return nodecover

        svfile = input[0]
        edgefile = input[1]
        combedge = parse_edge_file(edgefile)

        with open(svfile) as infile, open(output[0], "a") as outfile:

            for line in infile:
                # biallelic        1_165873        AltDel          497     2       s1,s2,s3        UCD,OBV         s1,s133016,s3   Angus,OBV
                linecomp = line.strip().split()
                sv_comp = linecomp[:3]
                ref_path = linecomp[5].split(",")
                nonref_path = linecomp[-2].split(",")
                ref_cover = extract_coverage(combedge, paths=ref_path)
                nonref_cover = extract_coverage(combedge, paths=nonref_path)

                print(*sv_comp, linecomp[5], linecomp[-2], ref_cover, nonref_cover, file=outfile)
