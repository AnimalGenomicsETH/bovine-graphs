#!/usr/bin/env python

rule repeat_mask:
    input:
        rules.combine_sv.output
    output:
        "analysis/bubble/{asb}_nonrefsv.fa.masked"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
            RepeatMasker -q -no_is -species cow {input}
        """


rule blastx_nonref:
    input: "analysis/bubble/{asb}_nonrefsv.fa.masked"
    output: "analysis/bubble/{asb}_nonrefsv_blastx.tsv"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "01:00"
    params:
        refprot = config["refprot"]
    shell:
        """

        diamond blastx --more-sensitive --db {params.refprot} --query {input} --threads {threads} \
                --out {output}.temp --evalue 1e-10 --max-target-seqs 1 --outfmt 6 \
                qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen salltitles

        # filter blastx result

        sort -k1,1 -k11,11nr -k12,12n  {output}.temp |
        sort -u -k1,1 --merge |
        awk '$4/$14>=0.7 && $3>=80' > {output} && rm {output}.temp

        """


def get_ref(assemb):
    refgenome = datstat.loc[datstat.assemb == assemb, "ascomp"].iloc[0].split(",")[0]
    return f"assembly/{refgenome}_full.fa"


rule create_extended_ref:
    input:
        lambda wildcards: get_ref(wildcards.asb),
        rules.repeat_mask.output
    output:
        "rna_seq/reference/{ref}+{asb}.fa"
    shell:
        """

        cat {input} > {output}

        """

rule generate_hisat_index:
    input:
        rules.create_extended_ref.output
    output:
        touch("rna_seq/reference/{asb}_{ref}_hisat_index_finished")
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """
            hisat2-build -p {threads} {input} rna_seq/reference/{wildcards.ref}+{wildcards.asb}
        """

rule map_transcriptome:
    input:
        rna1 = config["rna_basedir"] + "/{rna_anims}_qc_R1.fq.gz",
        rna2 = config["rna_basedir"] + "/{rna_anims}_qc_R2.fq.gz",
        refind = rules.generate_hisat_index.output
    output:
        "rna_seq/aligned/{ref}_{asb}/{rna_anims}_{asb}.bam"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00",
        disk_scratch = 50
    shell:
        """

        hisat2  --dta -x rna_seq/reference/{wildcards.ref}+{wildcards.asb} -1 {input.rna1} -2 {input.rna2} |
        samtools view -hu |
        samtools sort -T $TMPDIR -@ 10 -O BAM -o {output} -


        """

rule generate_hisat_linear:
    input:
        f"assembly/{config['reference']}_full.fa"
    output:
        touch(f"rna_seq/reference/{config['reference']}_hisat_index_finished")
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    params:
        reference = config["reference"]
    shell:
        """
            hisat2-build -p {threads} {input} rna_seq/reference/{params.reference}

        """

rule map_linear_transcriptome:
    input:
        rna1 = config["rna_basedir"] + "/{rna_anims}_qc_R1.fq.gz",
        rna2 = config["rna_basedir"] + "/{rna_anims}_qc_R2.fq.gz",
        refind = rules.generate_hisat_linear.output
    output:
        f"rna_seq/aligned/{config['reference']}/{{rna_anims}}_{config['reference']}.bam"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00",
        disk_scratch = 50
    params:
        reference = config["reference"]
    shell:
        """

        hisat2 --dta  -x rna_seq/reference/{params.reference} -1 {input.rna1} -2 {input.rna2} |
        samtools view -hu |
        samtools sort -T $TMPDIR -@ 10 -O BAM -o {output} -


        """

rule predict_gene_model:
    input:
        rules.repeat_mask.output
    output:
        "rna_seq/gene_model/{asb}_nonref_augustus.out",
        "rna_seq/gene_model/{asb}_nonref_augustus.gff"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    envmodules:
        "gcc/4.8.2", "gdc", "boost/1.55.0", "sqlite3/3.11", "gsl/1.16", "perl/5.18.4",
        "bamtools/2.5.1", "suitesparse/4.5.1", "openblas/0.2.13_seq", "augustus/3.3.2"
    shell:
        """


        augustus --species=human --UTR=on  {input} > {output[0]}

        grep -v "#" {output[0]} > {output[1]}

        """

rule analyze_gene_model:
    input:
        rules.predict_gene_model.output[1]
    output:
        "rna_seq/gene_model/{asb}_predict_summary.tsv"
    script: "../scripts/gene_model_analysis.R"


localrules: create_protein_fa
rule create_protein_fa:
    input: "rna_seq/gene_model/{asb}_nonref_augustus.out"
    output: "rna_seq/gene_model/{asb}_nonref_augustus_prot.fa"
    run:
        import re

        with open(input[0]) as infile, open(output[0], "a") as outfile:
            for line in infile:
                if line.startswith("# start gene"):
                    gene_id = line.strip().split()[-1]
                    contig_id = next(infile).strip().split()[0]
                    seqname = f">{gene_id}_{contig_id}"
                elif line.startswith("# protein sequence"):
                    if re.search(r"([A-Z]+)", line):
                        protseq = re.search(r"([A-Z]+)", line)[0]
                        nextline = next(infile)
                        while not nextline.startswith("# end gene"):
                            protseq += re.search(r"([A-Z]+)", nextline)[0]
                            nextline = next(infile)
                        else:
                            outfile.write(f"{seqname}\n{protseq}\n")


rule blastp_nonref:
    input: "rna_seq/gene_model/{asb}_nonref_augustus_prot.fa"
    output: "rna_seq/gene_model/{asb}_nonref_agustus_blastp.tsv"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "01:00"
    params:
        refprot = config["refprot"]
    shell:
        """

        #create protein fa from augustus results


        diamond blastp --more-sensitive --db {params.refprot} --query {input} --threads {threads} \
                --out {output}.temp --evalue 1e-10 --outfmt 6 \
                qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen salltitles

        # filter blastx result

        sort -k1,1 -k11,11nr -k12,12n  {output}.temp |
        sort -u -k1,1 --merge |
        awk '$4/$14>=0.7 && $3>=70' > {output}

        """

rule assemble_transcript:
    input:
        bam = "rna_seq/aligned/{ref}_{asb}/{rna_anims}_{asb}.bam",
        gffpred = "rna_seq/gene_model/{asb}_nonref_augustus.gff"
    output:
        "rna_seq/transcript_assembly/{asb}/{rna_anims}/{ref}+{asb}_{rna_anims}_temp_transcript",
        "rna_seq/transcript_assembly/{asb}/{rna_anims}/{ref}+{asb}_{rna_anims}_temp_abundance"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    params:
        gffinput = config["gffinput"]
    shell:
        """

        stringtie {input.bam} -G {input.gffpred} -G {params.gffinput} \
        -l {wildcards.asb}_{wildcards.rna_anims} -o {output[0]} -p {threads} -A {output[1]} -B

        """

rule merge_annotation:
    input:
        allanims = expand(
            "rna_seq/transcript_assembly/{{asb}}/{rna_anims}/{{ref}}+{{asb}}_{rna_anims}_temp_transcript", rna_anims=rna_anims),
        gffpred = "rna_seq/gene_model/{asb}_nonref_augustus.gff"
    output:
        "rna_seq/transcript_assembly/{asb}/{ref}+{asb}_anims.tsv",
        "rna_seq/transcript_assembly/{asb}/{ref}+{asb}_comb.gff"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    params:
        gffinput = config["gffinput"]
    shell:
        """

        echo {input.allanims} | tr ' ' '\\n' > {output[0]}

        stringtie --merge -G {input.gffpred} -G {params.gffinput} -l {wildcards.asb}_comb -o {output[1]} {output[0]}

        """

rule calculate_expression:
    input:
        bam = "rna_seq/aligned/{ref}_{asb}/{rna_anims}_{asb}.bam",
        gffmerged = "rna_seq/transcript_assembly/{asb}/{ref}+{asb}_comb.gff"
    output:
        "rna_seq/transcript_assembly/{asb}/{rna_anims}/{ref}+{asb}_{rna_anims}_merged_transcript",
        "rna_seq/transcript_assembly/{asb}/{rna_anims}/{ref}+{asb}_{rna_anims}_merged_abundance"
    threads: 10
    resources:
        mem_mb = 5000,
        walltime = "04:00"
    shell:
        """

        stringtie {input.bam} -e -G {input.gffmerged}  -o {output[0]} -p {threads} -A {output[1]} -B

        """

rule merge_expression:
    input:
        expand(
            "rna_seq/transcript_assembly/{{asb}}/{rna_anims}/{{ref}}+{{asb}}_{rna_anims}_merged_abundance", rna_anims=rna_anims)
    output: "rna_seq/transcript_assembly/{asb}/{ref}+{asb}_expression.tsv"
    run:

        from os.path import basename

        allfile = iter(input)

        def open_expression_file(procfile):
            anims = basename(procfile).split("_")[1]
            datin = pd.read_csv(procfile,
                                sep="\t",
                                skiprows=1,
                                header=None,
                                usecols=[0, 2, 4, 5, 8],
                                names=["gene_id", "contigs", "start_pos", "stop_pos", anims])
            # retain only on the non-reference contigs
            datin = datin[datin["contigs"].str.contains(r"^[m|b]\d+_\d+_\d+")]
            return datin

        # process first file
        combfile = open_expression_file(next(allfile))

        # for the subsequent file merge it with the first file

        for procfile in allfile:
            anims = basename(procfile).split("_")[1]
            infile = open_expression_file(procfile)[['gene_id', anims]]
            combfile = combfile.merge(infile, on="gene_id")

        # save the merged expression
        combfile.to_csv(str(output),
                        sep="\t",
                        index=False)
