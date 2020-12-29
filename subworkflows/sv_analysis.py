rule identify_bubble:
    input:
        "graph/{asb}_graph.gfa"
    output:
        "analysis/bubble/{asb}_bubble.tsv",
        "analysis/bubble/{asb}_biallelic_bubble.tsv",
        "analysis/bubble/{asb}_multiallelic_bubble.tsv",
        "analysis/bubble/{asb}_bubble.bed"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

        gfatools bubble {input} > {output[0]}

        awk '$5==2 {{ print $1,$2,$4,$5,$12 }}' {output[0]} > {output[1]}

        awk '$5>2 && $5 < 8 {{ print $1,$2,$4,$5,$12 }}' {output[0]} > {output[2]}

        awk '{{ print $1,$2,$2+1,$1"_"$2 }}' OFS="\t" {output[0]} > {output[3]}
        """

rule collect_biallelic_sv:
    input:
        "graph/{asb}_graph_len.tsv",
        "analysis/bubble/{asb}_biallelic_bubble.tsv"
    output:
        "analysis/bubble/{asb}_biallelic_sv.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """

            {workflow.basedir}/scripts/get_bialsv.py -a {wildcards.asb} > {output}

        """

rule extract_bialsv:
    input:
        "graph/{asb}_graph.gfa",
        rules.collect_biallelic_sv.output
    output:
        "analysis/bubble/{asb}_bialsv_seq.fa",
        "analysis/bubble/{asb}_bialsv_stat.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "00:30"
    shell:
        """

          {workflow.basedir}/scripts/get_bialseq.py -a {wildcards.asb}

        """

rule collect_multiallelic_sv:
    input:
        "graph/{asb}_graph_len.tsv",
        "analysis/bubble/{asb}_multiallelic_bubble.tsv",
        "analysis/colour_node/{asb}_nodecol.tsv"
    output:
        "analysis/bubble/{asb}_multiallelic_sv.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "01:00"
    shell:
        """
            {workflow.basedir}/scripts/get_multisv.py -a {wildcards.asb} > {output}
        """

rule extract_multisv:
    input:
        "graph/{asb}_graph.gfa",
        rules.collect_multiallelic_sv.output
    output:
        "analysis/bubble/{asb}_multisv_seq.fa",
        "analysis/bubble/{asb}_multisv_stat.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "00:30"
    shell:
        """
            {workflow.basedir}/scripts/get_multiseq.py -a {wildcards.asb}
        """

rule trace_paths:
    input:
        "remap/{asb}_edge_use.tsv",
        "analysis/bubble/{asb}_biallelic_sv.tsv",
        "analysis/bubble/{asb}_multiallelic_sv.tsv"
    output:
        "analysis/bubble/{asb}_path_trace.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    shell:
        """

        {workflow.basedir}/scripts/trace_path.py -g {wildcards.asb} -a {params.assemb} > {output}

        """

rule visualize_sv:
    input:
        rules.collect_biallelic_sv.output,
        rules.collect_multiallelic_sv.output,
        rules.construct_graph.output
    output: "analysis/bubble/{asb}_{svtype}_sv_viz.pdf"
    threads: 10
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    envmodules:
        "gcc/4.8.5",
        "graphviz/2.40.1",
        "python_cpu/3.7.4"
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    shell:
        """

            {workflow.basedir}/visualize/sv_viz.py -g {wildcards.asb} -c {params.assemb} -m {wildcards.svtype}

        """


rule combine_sv:
    input:
        "analysis/bubble/{asb}_bialsv_seq.fa",
        "analysis/bubble/{asb}_multisv_seq.fa"
    output:
        "analysis/bubble/{asb}_nonrefsv.fa"
    shell:
        """
            cat {input} > {output}
        """

rule combine_sv_woflank:
    input:
        bialfile = "analysis/bubble/{asb}_bialsv_stat.tsv",
        multifile = "analysis/bubble/{asb}_multisv_stat.tsv",
        graphfile = "graph/{asb}_graph.gfa"
    output:
        "analysis/bubble/{asb}_nonrefsv_woflank.fa"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        {workflow.basedir}/scripts/get_sv_woflank.py \
                -b {input.bialfile} \
                -m {input.multifile} \
                -g {input.graphfile} \
                -o {output}

        """


localrules: create_full_ref
rule create_full_ref:
    input:
        genome = f"assembly/{config['reference']}_full.fa",
        sv = "analysis/bubble/{asb}_nonrefsv.fa"
    output: f"extended_ref/{config['reference']}+{{asb}}.fa"
    shell:
        """

        cat {input.genome} {input.sv} > {output}

        """
localrules: create_breakpoint_bed
rule create_breakpoint_bed:
    input:
        bubble_file = "analysis/bubble/{asb}_bubble.tsv",
        bialsv_file = "analysis/bubble/{asb}_biallelic_sv.tsv",
        multisv_file = "analysis/bubble/{asb}_multiallelic_sv.tsv"
    output:
        left_bed = "analysis/bubble/{asb}_left_breakpoints.bed",
        right_bed = "analysis/bubble/{asb}_right_breakpoints.bed"
    run:
        import re

        # get right breakpoint
        right_bp = {}
        with open(input["bubble_file"]) as infile:
            for line in infile:
                chromo, left_side, right_side, *_ = line.strip().split()
                # only get the numeric part of the chromo
                chromo = [int(x) for x in chromo.split("_") if re.search(r"\d+", x)][0]
                right_bp[f"{chromo}_{left_side}"] = right_side

        def wrote_sv_bed(line, left_file, right_file, mutype="biallelic"):
            if mutype == "biallelic":
                # 1 165873 AltDel 497 2 s1 s2 s133016 s3
                line_comp = line.strip().split()
                chromo, leftcoord = line_comp[:2]
                start_node, stop_node = line_comp[5], line_comp[-1]
                sv_comp = f"{chromo}_{leftcoord}"
                leftcoord = int(leftcoord)
                # sv_comp, *sv_rest = line.strip().split()
                # _, *chromo, leftcoord = sv_comp.split("_")
                # only get the numeric part as the chromosome id
                # chromo = [int(x) for x in chromo if re.search(r"\d+", x)][0]
                # start_node, *_, stop_node = sv_rest[-2].split(",")
            if mutype == "multiallelic":
                # 1_535561        1898    5138    AltIns  s33,s60173,s36
                line_comp = line.strip().split()
                chromo, leftcoord = line_comp[0].split("_")
                leftcoord = int(leftcoord)
                sv_comp = f"{chromo}_{leftcoord}"
                start_node, *_, stop_node = line_comp[-1].split(",")
            # write the bed file
            left_file.write(
                (f"{chromo}\t{leftcoord-1}\t{leftcoord+1}\t{start_node}\t{stop_node}\t{sv_comp}\n"))
            svid = f"{chromo}_{leftcoord}"
            rightcoord = int(right_bp[svid])
            right_file.write(
                (f"{chromo}\t{rightcoord-1}\t{rightcoord+1}\t{start_node}\t{stop_node}\t{sv_comp}\n"))

        with open(input["bialsv_file"]) as bialfile, open(input["multisv_file"]) as multifile:
            with open(output["left_bed"], "a") as left_file, open(output["right_bed"], "a") as right_file:
                for line in bialfile:
                    # process the biallelic breakpoints
                    wrote_sv_bed(line, left_file, right_file, mutype="biallelic")
            # process the multiallelic breakpoints
                sv_processed = []
                mutlist = []
                for line in multifile:
                    svid = line.strip().split()[0]
                    # sv_comp, *sv_rest = line.strip().split()
                    # _, *chromo, leftcoord = sv_comp.split("_")
                    # # only get the numeric part as the chromosome id
                    # chromo = [int(x) for x in chromo if re.search(r"\d+", x)][0]
                    # svid = f"{chromo}_{leftcoord}"
                    if svid not in sv_processed:
                        sv_processed.append(svid)
                        wrote_sv_bed(line, left_file, right_file, mutype="multiallelic")

localrules: annotate_breakpoint
rule annotate_breakpoint:
    input:
        left_bed = rules.create_breakpoint_bed.output[0],
        right_bed = rules.create_breakpoint_bed.output[1]
    output:
        "analysis/bubble/{asb}_breakpoint_annot.tsv"
    params:
        gffinput = config["gffinput"]
    script: "../scripts/annot_breakpoints.py"


rule annot_sv:
    input:
        "analysis/bubble/{asb}_bubble.bed"
    output:
        "analysis/bubble/{asb}_bubble_annot.tsv"
    threads: 10
    resources:
        mem_mb = 1000,
        walltime = "10:00"
    params:
        gffinput = config["gffinput"]
    script:
        "../scripts/annot_sv.py"


rule visualize_exon:
    input:
        bialsv_file = "analysis/bubble/{asb}_biallelic_sv.tsv",
        multisv_file = "analysis/bubble/{asb}_multiallelic_sv.tsv",
        annot_file = "analysis/bubble/{asb}_breakpoint_annot.tsv"
    output: "analysis/bubble/{asb}_exon_viz.pdf"
    threads: 10
    params:
        assemb = lambda wildcards: get_assemb(wildcards.asb)
    envmodules:
        "gcc/4.8.5",
        "graphviz/2.40.1",
        "python_cpu/3.7.4"
    resources:
        mem_mb = 1000,
        walltime = "01:00"
    script: "../visualize/sv_viz_exon.py"
