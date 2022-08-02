#!/usr/bin/env python3
"""Prepare pc_translations new file, whitout RT, PAR_Y or redundancies

It needs gencode.v{version}.pc_translations.fa.gz and gencode.v{version}.annotation.gtf
"""

from argparse import ArgumentParser
import pandas as pd
import gzip
import os
import re
import subprocess


def convert_fasta_to_df(pc_transl):

    df_lines = []

    with gzip.open(pc_transl, 'rt') as f:

        seq = ""

        for line in f.readlines():

            line = line.strip()

            if line.startswith('>'):
                if seq != "":
                    header.append(seq)
                    df_lines.append(header)
                    seq = ""
                header = line[1:].split("|")
                header = [header[2], header[1], header[6], header[5]]

            else:
                line = line.replace("*", "X")
                seq = seq + line.strip()
        #Last line
        header.append(seq)
        df_lines.append(header)

    pc_trans_df = pd.DataFrame(df_lines,
                      columns = ['gene_ID', 'transcript_ID', 'gene_name', 'trasncript_name', 'sequence'])

    # Remove "_PAR_Y " transcripts
    pc_trans_df = pc_trans_df[pc_trans_df["transcript_ID"].str.contains("PAR_Y")==False]

    return pc_trans_df


def remove_RT_plague(annot_path, annot_prefix):
    
    gtf_header = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    annot_df = pd.read_csv(annot_path+"/"+annot_prefix+"_RT.gtf", sep="\t", names = gtf_header)
    annot_df['attribute']=annot_df['attribute'].str.split(";")
    
    RTtrans = []

    RT_list = open(annot_path+"/"+annot_prefix+".RT.list.csv", "w")
    RT_list.write("gene_name\tgene_id\ttranscript_id\n")
    RT_genes = open(annot_path+"/"+annot_prefix+".RT.genes.csv", "w")
    RT_genes.write("gene_name\tgene_id\n")

    for i in annot_df.index:
        if annot_df["feature"][i] == 'transcript':
            ttype = re.findall('^.*"(.*)"$', annot_df["attribute"][i][4])[0]
            gtype = re.findall('^.*"(.*)"$', annot_df["attribute"][i][2])[0]

            if gtype == "protein_coding" or gtype == "polymorphic_pseudogene" or gtype == "protein_coding_LoF":
                if ttype == "protein_coding" or ttype == "nonsense_mediated_decay" or ttype == "protein_coding_LoF"\
                or ttype == "non_stop_decay" or re.compile(r'TR_.*_gene').search(ttype)\
                or re.compile(r'IG_.*_gene').search(ttype):

                    transcript_id = re.findall('^.*"(.*)"$', annot_df["attribute"][i][1])[0]
                    RTtrans.append(transcript_id)
                    gene = re.findall('^.*"(.*)"$', annot_df["attribute"][i][3])[0]
                    gene_id = re.findall('^.*"(.*)"$', annot_df["attribute"][i][0])[0]

                    RT_list.write(gene+"\t"+gene_id+"\t"+transcript_id+"\n")
                    RT_genes.write(gene+"\t"+gene_id+"\n")

    RT_list.close()
    RT_genes.close()

    pc_trans_df_noRT = pc_trans_df[~pc_trans_df['transcript_ID'].isin(RTtrans)].reset_index(drop=True)

    return pc_trans_df_noRT


def convert_df_to_fasta(pc_trans_df_noRT_u):
    line_len = 60

    with gzip.open(pc_transl_path+"/"+pc_transl_prefix+".NoRT.u.fa.gz", 'wt') as f:

        for i in pc_trans_df_noRT_u.index:
            header = "|".join([">X", pc_trans_df_noRT_u["transcript_ID"][i], pc_trans_df_noRT_u["gene_ID"][i], \
                               "-", "-", pc_trans_df_noRT_u["trasncript_name"][i], pc_trans_df_noRT_u["gene_name"][i], \
                               str(len(pc_trans_df_noRT_u["sequence"][i]))])
            seq = pc_trans_df_noRT_u["sequence"][i]
            chunks = [seq[i:i+line_len] for i in range(0, len(seq), line_len)]

            f.write(header+"\n")

            for j in chunks:
                f.write(j+"\n")
                
    return 1


if __name__ == "__main__":

    ap = ArgumentParser(description=__doc__)
    ap.add_argument("pc_transl_file",
                    help="input pc_translations .FA.GZ file")
    ap.add_argument("annot_file",
                    help="input annotations GTF.GZ file")
    args = ap.parse_args()

    pc_transl_file = args.pc_transl_file
    annot_file = args.annot_file



    annot_path = os.path.dirname(annot_file)
    annot_prefix = ".".join(os.path.basename(annot_file).split(".")[:-2])

    pc_transl_path = os.path.dirname(pc_transl_file)
    pc_transl_prefix = ".".join(os.path.basename(pc_transl_file).split(".")[:-2])


    ### Convert translations fasta to DF
    pc_trans_df = convert_fasta_to_df(pc_transl_file)


    ### Create GTF with RT
    search_command = ['zgrep', '"eadthro"', annot_file]
    subprocess.call(' '.join(search_command), stdout=open(annot_path+"/"+annot_prefix+"_RT.gtf", 'w'), shell=True)


    ### Remove RT from translations DF
    pc_trans_df_noRT = remove_RT_plague(annot_path, annot_prefix)


    ### Make it non-redundant
    pc_trans_df_noRT['gene_seq'] = pc_trans_df_noRT["gene_ID"]+"_"+pc_trans_df_noRT["sequence"]

    # DF with duplicates
    #df_duplicated = pc_trans_df_noRT[pc_trans_df_noRT.duplicated(subset="gene_seq", keep='first')]

    # DF_clean
    pc_trans_df_noRT_u = pc_trans_df_noRT.drop_duplicates(subset=['gene_seq']).reset_index(drop=True)


    ### Convert DF to fasta file again
    result = convert_df_to_fasta(pc_trans_df_noRT_u)
