#!/usr/bin/env python

import pandas as pd
import sys
configfile: "config/config.yaml"
workdir: config["workdir"]


datstat = pd.read_csv("config/graph_comp.tsv",
                      sep=" ",
                      header=None,
                      names=["assemb", "ascomp"])

graphcon = list(datstat["assemb"])
svlist = ["biallelic", "multiallelic"]
reflist = [x.split(",")[0] for x in datstat.loc[:, "ascomp"]]

# jobs without submission to cluster
localrules: combine_sv, create_extended_ref, analyze_gene_model, merge_expression

# parse output file
include: "subworkflows/pipeline_output.py"
core_out = core_analysis_output()
sv_out = sv_analysis_output()
rna_anims, rna_out = rna_analysis_output(include_rna_pipeline=config["rna_seq"])
dna_out = wgs_analysis_output()

rule all:
    input:
        expand("reports/{asb}_report.pdf", asb=graphcon)


def get_assemb(assemb):
    allcomp = datstat.loc[datstat.assemb == assemb, "ascomp"].iloc[0].split(",")
    return allcomp


rule construct_graph:
    input:
        lambda wildcards: [f"assembly/{x}.fa" for x in get_assemb(wildcards.asb)]
    output:
        "graph/{asb}_graph.gfa",
        "graph/{asb}_graph_len.tsv",
        "graph/{asb}_graph_link.tsv"
    threads: 10
    resources:
        mem_mb = 12000,
        walltime = "02:00"
    shell:
        """

        minigraph --inv no -xggs -t {threads} {input}  > {output[0]}

        awk '$1~/S/ {{ split($5,chr,":"); split($6,pos,":"); split($7,arr,":");
            print $2,length($3),chr[3],pos[3],arr[3] }}' {output[0]} > {output[1]}

        awk '$1 == "L"' {output[0]} > {output[2]}

        """

rule remap_graph:
    input:
        rules.construct_graph.output[0],
        "assembly/{anims}.fa"
    output:
        "remap/{asb}/{anims}_{asb}.gaf"
    threads: 10
    resources:
        mem_mb = 10000,
        walltime = "04:00"
    shell:
        """
         minigraph -t {threads} --cov -x asm {input[0]} {input[1]} > {output}

        """

rule comb_coverage:
    input:
        lambda wildcards: [f"remap/{wildcards.asb}/{x}_{wildcards.asb}.gaf" for x in get_assemb(wildcards.asb)]
    output:
        "remap/{asb}_coverage.tsv",
        "remap/{asb}_edge_use.tsv"
    params:
        anims = lambda wildcards: get_assemb(wildcards.asb)
    threads: 5
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

            {workflow.basedir}/scripts/comb_coverage.py -g {wildcards.asb} -a {params.anims}

        """

rule colour_node:
    input:
        rules.construct_graph.output[1],
        ancient(rules.comb_coverage.output[0])
    output:
        "analysis/colour_node/{asb}_nodecol.tsv",
        "analysis/colour_node/{asb}_nodemat.tsv"
    threads: 5,
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    shell:
        """
            {workflow.basedir}/scripts/colour_node.R {wildcards.asb} {params.assemb}
        """

localrules: identify_core_nonref
rule identify_core_nonref:
    input:
        rules.construct_graph.output[1],
        rules.colour_node.output[1]
    output:
        "analysis/core_nonref/{asb}_core_analysis.tsv",
        "analysis/core_nonref/{asb}_nonref_analysis.tsv",
        multiext("analysis/core_nonref/{asb}_core_flex_sim", ".tsv", ".png", ".pdf"),
        multiext("analysis/core_nonref/{asb}_nonref_shared_count", ".png", ".pdf"),
        multiext("analysis/core_nonref/{asb}_nonref_shared_len", ".png", ".pdf")
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    script:
        "scripts/run_core_nonref.R"

# Add workflow for genetic distance analysis
include: "subworkflows/mash_distance.py"

# Add workflow for sv analysis
include: "subworkflows/sv_analysis.py"

if config["rna_seq"]:
    # Add workflow for functional analysis
    include: "subworkflows/rnaseq_analysis.py"

# Add workflow for wgs analysis
if config["dna_seq"]:
    include: "subworkflows/wgs_analysis.py"

rule generate_report:
    input:
        core_out,
        sv_out,
        rna_out,
        dna_out
    output:
        "reports/{asb}_report.pdf"
    threads: 10
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb),
        include_rna = config["rna_seq"]
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """
            {workflow.basedir}/reports/generate_report.py -a {wildcards.asb} -r {params.include_rna} -c {params.assemb} 
        """
