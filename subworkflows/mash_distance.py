#!/usr/bin/env python

import glob
import re

assemblies = [re.findall(r"/([A-Za-z]+).fa", x)[0] for x in glob.glob("assembly/*fa") if "full" not in x]

rule sketch_assembly:
    input: "assembly/{assemb}.fa"
    output: "tree/{assemb}.fa.msh"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        mash sketch -p {threads} -o {output} {input}

        """

localrules: combined_sketch
rule combined_sketch:
    input: expand("tree/{assemb}.fa.msh", assemb=assemblies)
    output: "tree/combined_sketch.msh"
    shell:
        """

        mash paste {output} {input}

        """

rule estimate_distance:
    input: "tree/combined_sketch.msh"
    output: "tree/combined_distance.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

        mash dist {input} {input} > {output}

        """

rule visualize_tree:
    input: "tree/combined_distance.tsv"
    output:
        "tree/pairwise_distance.tsv",
        "tree/assembly_phylo_tree.pdf"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "10:00"
    params:
        ref = config["reference"]
    script: "../scripts/phylo_tree_assembly.R"
