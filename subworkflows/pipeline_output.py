#!/usr/bin/env python


def core_analysis_output():
    """

    Construct list for core genome analysis output

    """
    core_out = []

    core_out.extend(expand("analysis/colour_node/{asb}_nodecol.tsv", asb=graphcon))
    core_out.extend(expand("analysis/colour_node/{asb}_nodemat.tsv", asb=graphcon))
    core_out.extend(expand("analysis/core_nonref/{asb}_core_analysis.tsv", asb=graphcon))
    # phylogenetic tree across assemblies
    core_out.extend(["tree/assembly_phylo_tree.pdf"])

    return core_out


def sv_analysis_output():
    """

    Construct list for structural variations output

    """

    sv_out = []
    sv_out.extend(expand("analysis/bubble/{asb}_nonrefsv.fa", asb=graphcon))
    sv_out.extend(expand("analysis/bubble/{asb}_nonrefsv_woflank.fa", asb=graphcon))
    sv_out.extend(expand("analysis/bubble/{asb}_breakpoint_annot.tsv", asb=graphcon))
    sv_out.extend(expand("analysis/bubble/{asb}_{svtype}_sv_viz.pdf", asb=graphcon, svtype=svlist))
    sv_out.extend(expand("analysis/bubble/{asb}_exon_viz.pdf", asb=graphcon))

    # trace paths in the bubbles
    sv_out.extend(expand("analysis/bubble/{asb}_path_trace.tsv", asb=graphcon))

    # full extended reference
    sv_out.extend(expand(f"extended_ref/{config['reference']}+{{asb}}.fa", asb=graphcon))

    return sv_out


def rna_analysis_output(include_rna_pipeline=True):
    """
    Construct output for the rna_analysis pipeline

    Input: config["rna_seq"] if True will include all rna analysis output (default)
    Output: None if rna_seq not included otherwise all rna_seq analysis file output

    """

    rna_out = []
    rna_anims = []
    if include_rna_pipeline:

        rna_out.extend(expand("analysis/bubble/{asb}_nonrefsv.fa.masked", asb=graphcon))
        rna_out.extend(expand("rna_seq/reference/{ref}+{asb}.fa", zip, ref=reflist, asb=graphcon))

        # add wildcards from animal in transcriptome
        rna_anims, = glob_wildcards(f"{config['rna_basedir']}/{{rna_anims}}_qc_R1.fq.gz")
        if not rna_anims:
            sys.exit("No transcriptome data found. Possibly path is incorrect")

        rna_out.extend(expand(
            [f"rna_seq/transcript_assembly/{asb}/{{rna_anims}}/{ref}+{asb}_{{rna_anims}}_merged_transcript" for ref, asb in zip(reflist, graphcon)], rna_anims=rna_anims))

        # add result for the gene_prediction
        rna_out.extend(expand("rna_seq/gene_model/{asb}_predict_summary.tsv", asb=graphcon))

        # add mapping to reference results

        rna_out.extend(expand(
            f"rna_seq/aligned/{config['reference']}/{{rna_anims}}_{config['reference']}.bam", rna_anims=rna_anims))

        # add merge expression results

        rna_out.extend([f"rna_seq/transcript_assembly/{asb}/{ref}+{asb}_expression.tsv"
                        for ref, asb in zip(reflist, graphcon)])
        # add blast results

        rna_out.extend(expand("analysis/bubble/{asb}_nonrefsv_blastx.tsv", asb=graphcon))
        rna_out.extend(expand("rna_seq/gene_model/{asb}_nonref_agustus_blastp.tsv", asb=graphcon))

    return rna_anims, rna_out


def wgs_analysis_output(include_dna=config["dna_seq"]):

    dna_out = []

    if include_dna:
        dna_out.append("wgs/stat/wgs_mapping_stat.tsv")

    return dna_out
