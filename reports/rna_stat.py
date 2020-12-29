#!/usr/bin/env python


def gene_model_string(assembly):
    model_string = "| Gene model features | Statistics |\n"
    model_string += "|------------------ | ----- |\n"
    with open(f"rna_seq/gene_model/{assembly}_predict_summary.tsv") as infile:
        # skip header
        next(infile)
        for line in infile:
            feature, stat = line.strip().split(";")
            model_string += f"|{feature}|{stat}|\n"
    return model_string
