#!/usr/bin/env python

import csv


pep_cutoff = snakemake.params.pep_cutoff
seq_db = snakemake.wildcards.sdb


if seq_db == "gencode":
    prot_to_gene = dict()
    with open(snakemake.input.meta_file,"r") as mfh:
        reader = csv.DictReader(mfh,dialect="excel-tab")
        for counter,row in enumerate(reader,start=1):
            if row["db_name"] == "gencode":
                protter_id = "seq{}#gencode#".format(counter)
                gene_id = row["db_seq_id"].split("|")[2]
                bare_gene_id = gene_id.split(".",1)[0]
                prot_to_gene[protter_id] = bare_gene_id


with open(snakemake.input.psms,"r") as ifh:
    reader = csv.DictReader(ifh,dialect="excel-tab")
    out_headings = reader.fieldnames
    if seq_db == "gencode":
        out_headings += ["gene id"]
    with open(snakemake.output.psms,"w") as ofh:
        writer = csv.DictWriter(ofh,out_headings,dialect="excel-tab")
        writer.writeheader()
        for row in reader:
            pep_score = float(row["percolator PEP"])
            if pep_score <= pep_cutoff:
                if seq_db == "gencode":
                    prot_ids = row["protein id"].split(",")
                    gene_ids = [prot_to_gene[x] if x in prot_to_gene else "NA"
                                for x in prot_ids]
                    row["gene id"] = ",".join(gene_ids)
                writer.writerow(row)
