#!/usr/bin/env python3
"""Process output PSM files to get the proteo CSV file

It needs protter_psm_output.zip and the config.yaml file from the same running
"""

from argparse import ArgumentParser
from Bio import SeqIO
from collections import Counter, defaultdict
from zipfile import ZipFile
import glob
import gzip
import os
import pandas as pd
import shutil
import yaml

from common import load_config_file

###################
#def load_config_file(config_file):
#    with open(config_file) as f:
#        return yaml.safe_load(f)
##################################################### 

# Get dataset data
def get_dataset_data(psm_path, db, ds):

    db_ds_folder = psm_path+"/"+db+"/"+ds
    list_proc = sorted(glob.glob(db_ds_folder+"/*gencode*"))

    ds_df = pd.DataFrame(columns=['sequence', 'expt', 'percolator q-value', 'percolator PEP', \
                                  'gene id', 'tissue'])

    for proc in list_proc:

        expt = os.path.basename(proc).split(".")[3]
        expt_df = pd.read_csv(proc, sep="\t", header = 0)
        expt_df['sequence'] = expt_df['sequence'].str.replace(r'[0-9.\[\]]+', '', regex=True)
        expt_df['expt'] = expt

        if not 'gene id' in expt_df.columns:
            print("\tNo genes in ",proc)
            expt_df['gene id'] = ""
        if not 'tissue' in expt_df.columns:
            print("\tNo tissue in ",proc)
            expt_df['tissue'] = ""
        filter_exp_df = expt_df[['sequence', 'expt', 'percolator q-value', 'percolator PEP', \
                                 'gene id', 'tissue']]

        ds_df = pd.concat([ds_df, filter_exp_df])

    ds_df = ds_df.reset_index(drop = True)
    return ds_df
# END Get dataset data


# Make Unique. Counting experiments and tissues
def make_unique(in_tsv, out_tsv):
    with open(in_tsv, 'r') as fp, open(out_tsv, 'w') as out:
        # Discard header
        fp.readline()

        # Add header to new file
        out.write("sequence\tpercolator PEP\tgene id\tcount\ttag\ttissues\n")

        # Inicialice variables with first line
        line1 = fp.readline().strip().split("\t")  

        old_expt = line1[1]
        old_pep = line1[0]
        count = 1
        best = float(line1[3])
        old_genes = line1[4]
        tag = ""
        if not (len(set(old_genes))==1):
                tag = "Moonlighting"
        tissues = [line1[5]]

        for line in fp:
            line_split = line.strip().split("\t")

            if old_pep != line_split[0]:

                tissues_count = Counter(tissues)
                tissues_dict = {}
                for key in tissues_count:  
                    val=tissues_count[key]
                    tissues_dict[key] = val

                tag = ""
                genes = old_genes.split(",")
                if not (len(set(genes))==1):
                    tag = "Moonlighting"

                tissues_str = str(tissues_dict).replace(',', ';')
                out.write(old_pep+"\t"+str(best)+"\t"+genes[0]+"\t"+str(count)+"\t"+tag+"\t"+tissues_str+"\n")

                count = 0
                best = 1
                old_expt = 0
                tissues = []

            if line_split[1] != old_expt:
                count += 1
                tissues.append(line_split[5])

            if float(line_split[3]) < best:
                best = float(line_split[3])

            old_pep = line_split[0]
            old_expt = line_split[1]
            old_genes =line_split[4]

        tissues_count = Counter(tissues)
        tissues_dict = {}
        for key in tissues_count:  
            val=tissues_count[key]
            tissues_dict[key] = val

        tag = ""
        genes = old_genes.split(",")
        if not (len(set(genes))==1):
            tag = "Moonlighting"

        tissues_str = str(tissues_dict).replace(',', ';')
        out.write(old_pep+"\t"+str(best)+"\t"+genes[0]+"\t"+str(count)+"\t"+tag+"\t"+tissues_str+"\n")
# END Make Unique


# Extract peps hash
def extract_peps_hash(transl_file, in_tsv, out_tsv):
    geneid = {}
    genelist = {}
    protid = {}
    transcripts = []
    old_gene = ""

    with gzip.open(transl_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    # Inicialice variables with first line
    line1 = records[0].id.split("|")
    old_gene = line1[2].split('.')[0]
    transcripts.append(line1[1])
    geneid[line1[1]] = old_gene
    prot = str(records[0].seq).replace('L', 'I')
    protid[line1[1]] = prot

    for i in range(1,len(records)):
        chunks = records[i].id.split("|")
        gene = chunks[2].split('.')[0]

        if old_gene != gene:
            genelist[old_gene] = transcripts
            transcripts = []

        transcripts.append(chunks[1])
        geneid[chunks[1]] = gene
        prot = str(records[i].seq).replace('L', 'I')
        protid[chunks[1]] = prot
        old_gene = gene

    genelist[old_gene] = transcripts


    with open(in_tsv, 'r') as in_file, open(out_tsv, 'w') as out_file:

        #Remove header
        in_file.readline()

        # Insert headerin new file
        out_file.write("peptide\tgeneid\tmatched\tnot matched\tpercolator PEP\tcount\ttissues\n")

        for line in in_file:
            chunks = line.strip().split("\t")
            peptide = chunks[0]
            haspep = []
            notpep = []
            geneid = chunks[2]
            if chunks[4] == "Moonlighting":
                continue
            trans = genelist[geneid]
            for i in range(len(trans)):
                if peptide in protid[trans[i]]:
                    haspep.append(trans[i])
                else:
                    notpep.append(trans[i])
            haspep = ";".join(haspep)
            notpep = ";".join(notpep)

            if haspep == '':
                out_file.write(peptide+"\tNot in gencode\t\t\t"+chunks[1]+"\t"+chunks[3]+"\t"+chunks[5]+"\n")
            else:
                out_file.write(peptide+"\t"+geneid+"\t"+haspep+"\t"+notpep+"\t"+chunks[1]+"\t"+chunks[3]+"\t"+chunks[5]+"\n")
# END Extract peps hash


# Filter semitryptics
def filter_semitryptics(transl_file, in_tsv, out_tsv):
    transseq = {}

    with gzip.open(transl_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):

            chunks = record.id.split("|")
            enstrans = chunks[1]
            ensgene = chunks[2]
            seq = str(record.seq).replace('L', 'I')
            transseq[enstrans] = seq

    with open(in_tsv, 'r') as in_file, open(out_tsv, 'w') as out_file:

        #Remove header
        in_file.readline()

        # Insert headerin new file
        out_file.write("peptide\tgeneid\tmatched\tnot matched\tpercolator PEP\tcount\ttissues\n")

        for line in in_file:
            chunks = line.strip().split("\t")
            peptide = chunks[0].replace('L', 'I')
            haspep = chunks[2].split(";")
            typeN = "Semi"
            letters = [*peptide]

            # Remove semi-tryps
            if not (letters[-1] == "R" or letters[-1] == "K"):
                continue

            for i in haspep:
                transseq[i] = transseq[i].replace('L','I')
                letters = transseq[i].split(peptide)
                if letters[0] == "X":
                    continue

            # Remove peptides with 1 PSM, wrong PEP score or wrong lenght
            peplen = len(peptide)
            PEPscore = float(chunks[4])
            psmcount = int(chunks[5])
            if psmcount < 2 or PEPscore > 0.001 or peplen > 35 or peplen < 7:
                continue

            # Remove peptides with more than 1 missed cleavage
            letters = [*peptide]
            count = 0

            for i in range(0,len(letters)-1):
                if (letters[i] == "K" and letters[i+1] != "P") or (letters[i] == "R" and letters[i+1] != "P"):
                    count += 1
            if count > 1:
                continue

            out_file.write(line)
# END Filter semitryptics

            
### MAIN PROGRAM

if __name__ == "__main__":

    ap = ArgumentParser(description=__doc__)
    ap.add_argument("config_file",
                    help="Input config YAML file indicating the database(s) and"
                         " dataset(s) for which PSM files have been gathered.")
    ap.add_argument("psm_file",
                    help="Input ZIP file containing the gathered PSM files.")
    args = ap.parse_args()

    config_file = args.config_file
    psm_file = args.psm_file

    config = load_config_file(config_file)
    psm_path = os.path.dirname(os.path.abspath(psm_file))

    ### Get database(s) and datasets
    db_to_datasets = defaultdict(list)
    for ds in config["datasets"]:
        if not config["datasets"][ds].get("enabled",False):
            continue
        for db in config["datasets"][ds]["dbs"]:
            if not config["dbs"][db].get("enabled",False):
                continue
            db_to_datasets[db].append(ds)

    ### Unzip gather_psm_files output        
    with ZipFile(psm_file, 'r') as psm_zip:
        # Extract all the contents of zip file in current directory
        print("Extracting file "+psm_file+"...")
        psm_zip.extractall(path=psm_path)


    for db,datasets in db_to_datasets.items():
        print("{} RUNNING...".format(db))
        
        ### List all PSMs for each experiment and sort (by peptide, k1) all results together in one big DF
        all_df = pd.DataFrame(columns=['sequence', 'expt', 'percolator q-value', 'percolator PEP', \
                                    'gene id', 'tissue'])
        for ds in datasets:
            print("\tExtract dataset '{}' from '{}'...".format(ds,db))       
            all_df = pd.concat([all_df, get_dataset_data(psm_path, db, ds)], ignore_index=True)
            
        # Sorting PSMs from DB
        print("\tSorting PSMs for '{}'...".format(db))  
        all_sort_df = all_df.sort_values(by=['sequence', 'expt'], ignore_index = True)
        all_sort_df.to_csv(os.path.join(psm_path, db, 'All.sort.tsv'), sep="\t", header=True, index=False)

        # Making unique
        print("\tMaking PSMs for '{}' unique...".format(db))  
        make_unique(os.path.join(psm_path, db, 'All.sort.tsv'), os.path.join(psm_path, db, 'All.sortu.tsv'))

        transl_file = config['dbs'][db]['paths']['gencode']

        # Extracting peps
        print("\tExtracting pep matches for '{}'...".format(db))  
        extract_peps_hash(transl_file, os.path.join(psm_path, db, 'All.sortu.tsv'), os.path.join(psm_path, db, 'All.sortu.out2.tsv'))

        # Filtering peptides
        print("\tFiltering pepides for '{}'...".format(db))  
        filter_semitryptics(transl_file, os.path.join(psm_path, db, 'All.sortu.out2.tsv'), os.path.join(psm_path, db, 'All.semitryps.tsv'))

        # Sort by gene id, save final proteo file and remove intermediate files
        print("\tCreating PROTEO file for '{}'...".format(db))  
        allsemitryps = pd.read_csv(os.path.join(psm_path, db, 'All.semitryps.tsv'), sep = "\t", header = 0)
        proteo_file = allsemitryps.sort_values('geneid').reset_index(drop=True)
        proteo_file.to_csv(os.path.join(psm_path,'proteo_'+db+'.csv'), sep = ",", header = True, index = False)
        print("{} FINISHED\n".format(db))
        
        # Remove unzipped protter_psm folder
        try:
            shutil.rmtree(os.path.join(psm_path,db))
        except:
            print("Error: {} folder could not be removed".format(db))
