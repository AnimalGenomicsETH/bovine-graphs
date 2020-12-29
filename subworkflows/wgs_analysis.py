#!/usr/bin/env python

fastqdir = config["dna_basedir"]


def get_dna_ref(combref=config["dna_ref"]):
    refmap = {}
    for ref in combref:
        if ref in graphcon:
            refmap[f"{ref}pan"] = f"extended_ref/{config['reference']}+{ref}.fa"
        else:
            refmap[f"{ref}lin"] = f"assembly/{ref}_full.fa"
    return refmap


refmap = get_dna_ref()

localrules: link_reference
rule link_reference:
    input: lambda wildcards: refmap[wildcards.dnaref]
    output: "wgs/reference/{dnaref}.fa"
    shell:
        """

        ln -s $PWD/{input} $PWD/{output}

        """

rule index_wgs:
    input: "wgs/reference/{dnaref}.fa"
    output: touch("wgs/reference/{dnaref}_index_finished")
    threads: 18
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    envmodules:
        "gcc/4.8.5",
        "bwa"
    shell:
        """
            bwa index {input}
        """


rule map_genome:
    input:
        r1 = fastqdir + "/{dna_anims}_R1.fastq.gz",
        r2 = fastqdir + "/{dna_anims}_R2.fastq.gz",
        ref = "wgs/reference/{dnaref}.fa",
        refindex = rules.index_wgs.output
    output:
        bam = "wgs/bam/{dna_anims}_{dnaref}.bam"
    shadow: "full"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "24:00",
        disk_scratch = 50
    params:
        anims = "{dna_anims}",
        assemb = "{dnaref}",
        fastqdir = fastqdir
    envmodules:
        "gcc/4.8.5",
        "bwa"
    shell:
        """

        fastp --thread {threads} -i {input.r1} -I {input.r2} -o R1_qc.fq.gz -O R2_qc.fq.gz

        bwa mem -t {threads} -R "@RG\\tID:anim1\\tPL:Illumina\\tLB:lib1\\tSM:{params.anims}" {input.ref} R1_qc.fq.gz R2_qc.fq.gz |
        samblaster|
        samtools sort -@ {threads} -T $TMPDIR -O BAM -o {output} -

        """

rule index_bam:
    input: "wgs/bam/{dna_anims}_{dnaref}.bam"
    output: "wgs/bam/{dna_anims}_{dnaref}.bam.bai"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    shell:
        """

        samtools index -@ {threads} {input} 

        """

rule stat_mapping:
    input:
        bam = rules.map_genome.output,
        bai = rules.index_bam.output
    output:
        "wgs/stat/{dnaref}/{dna_anims}_{dnaref}_stat.tsv"
    threads: 10
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    params:
        anims = "{dna_anims}",
        assemb = "{dnaref}"
    shell:
        """

        {workflow.basedir}/scripts/collect_bam_stat.py -a {params.anims} -m {params.assemb} > {output}

        """

localrules: collect_stat
rule collect_stat:
    input:
        expand("wgs/stat/{dnaref}/{dna_anims}_{dnaref}_stat.tsv",
               dna_anims=config["dna_anims"], dnaref=refmap.keys())
    output: "wgs/stat/wgs_mapping_stat.tsv"
    shell:
        """
        (echo -e 'totmap\\tunmap\\tmapq0\\tmapq10\\tmapq60\\tal99\\talp\\tclip';cat {input}) > {output}

        """
