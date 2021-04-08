import csv

from common import get_group_enzyme,get_group_meta_value,load_sample_sheet


input_psm_file = snakemake.input.psms
output_psm_file = snakemake.output.psms
target_meta_file = snakemake.input.target_meta_file
sample_meta_file = snakemake.input.sample_meta_file

pep_cutoff = snakemake.params.pep_cutoff
seq_db = snakemake.wildcards.sdb
ds = snakemake.wildcards.ds
subset = snakemake.wildcards.subset
grouping = snakemake.wildcards.grouping
group = snakemake.wildcards.group

default_enzyme = snakemake.config["software"]["percolator"]["default_enzyme"]

samples = load_sample_sheet(sample_meta_file)
enzyme = get_group_enzyme(ds,subset,grouping,group,samples,default=default_enzyme)
tissue = get_group_meta_value(ds,subset,grouping,group,samples,"tissue")


if seq_db in ("gencode","upstream"):
    prot_to_gene = dict()
    with open(target_meta_file,"r") as mfh:
        reader = csv.DictReader(mfh,dialect="excel-tab")
        for counter,row in enumerate(reader,start=1):
            if seq_db == "gencode" and row["db_name"] == seq_db:
                protter_id = "seq{}#gencode#".format(counter)
                db_seq_id = row["db_seq_desc"].split(maxsplit=1)[0]
                gene_id = db_seq_id.split("|")[2]
                bare_gene_id = gene_id.split(".",1)[0]
                prot_to_gene[protter_id] = bare_gene_id
            elif seq_db == "upstream" and row["db_name"] == seq_db:
                protter_id = "seq{}#upstream#".format(counter)
                db_seq_id = row["db_seq_desc"].split(maxsplit=1)[0]
                gene_id = db_seq_id.split("|")[0]
                bare_gene_id = gene_id.split(".",1)[0]
                prot_to_gene[protter_id] = bare_gene_id


with open(input_psm_file,"r") as ifh:
    reader = csv.DictReader(ifh,dialect="excel-tab")
    out_headings = reader.fieldnames + ["enzyme"]
    if seq_db in ("gencode","upstream"):
        out_headings.append("gene id")
    if tissue is not None:
        out_headings.append("tissue")
    with open(output_psm_file,"w") as ofh:
        writer = csv.DictWriter(ofh,out_headings,dialect="excel-tab")
        writer.writeheader()
        for row in reader:
            pep_score = float(row["percolator PEP"])
            if pep_score <= pep_cutoff:
                row["enzyme"] = enzyme
                if seq_db in ("gencode","upstream"):
                    prot_ids = row["protein id"].split(",")
                    gene_ids = [prot_to_gene[x] if x in prot_to_gene else "NA"
                                for x in prot_ids]
                    row["gene id"] = ",".join(gene_ids)
                if tissue is not None:
                    row["tissue"] = tissue
                writer.writerow(row)
