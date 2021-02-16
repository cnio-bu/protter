import os
import re
import shutil
from tempfile import TemporaryDirectory

import pandas as pd
import requests


def enhance_sample_metadata(ds_tab):

    tissue_map = {
        ("Adult","Adrenal gland"): "adrenal gland",
        ("Adult","B cells"): "B cells",
        ("Adult","CD4 Tcells"): "CD4+ T cells",
        ("Adult","CD8 Tcells"): "CD8+ T cells",
        ("Adult","CD8  Tcells"): "CD8+ T cells",
        ("Adult","Colon"): "colon",
        ("Adult","Esophagus"): "esophagus",
        ("Adult","Frontal cortex"): "frontal cortex",
        ("Adult","Gallbladder"): "gallbladder",
        ("Adult","Heart"): "heart",
        ("Adult","Kidney"): "kidney",
        ("Adult","Liver"): "liver",
        ("Adult","Lung"): "lung",
        ("Adult","Monocytes"): "monocytes",
        ("Adult","NK cells"): "NK cells",
        ("Adult","Ovary"): "ovary",
        ("Adult","Pancreas"): "pancreas",
        ("Adult","Platelets"): "platelets",
        ("Adult","Prostate"): "prostate",
        ("Adult","Rectum"): "rectum",
        ("Adult","Retina"): "retina",
        ("Adult","Spinal cord"): "spinal cord",
        ("Adult","Testis"): "testis",
        ("Adult","Urinary bladder"): "urinary bladder",
        ("Fetus","Brain"): "fetal brain",
        ("Fetus","Gut"): "fetal gut",
        ("Fetus","Heart"): "fetal heart",
        ("Fetus","Liver"): "fetal liver",
        ("Fetus","Ovary"): "fetal ovary",
        ("Fetus","Placenta"): "fetal placenta",
        ("Fetus","Testis"): "fetal testis"
    }

    base_url = "http://ftp.pride.ebi.ac.uk/pride/data/archive"
    meta_file_name = "Metadata_Draft_map_of_human_proteoms_kim_et_al.xls"
    meta_file_url = "{}/2014/04/PXD000561/{}".format(base_url,meta_file_name)

    with TemporaryDirectory() as tmp_dir:

        local_meta_file = os.path.join(tmp_dir,meta_file_name)
        with requests.get(meta_file_url,stream=True) as r:
            with open(local_meta_file,"wb") as out_f:
                shutil.copyfileobj(r.raw,out_f)

        df = pd.read_excel(local_meta_file,sheet_name="Sheet1",
                           usecols=["Experiment name","Raw File name",
                                    "Developmental stage","Sample","Enzyme"],
                           engine="xlrd")
        df.dropna(inplace=True)  # drop copyright notice from data

    study_meta = dict()
    for _,row in df.iterrows():
        raw_file_name = row["Raw File name"]
        experiment = row["Experiment name"]
        sample_stage = row["Developmental stage"]
        sample_tissue = row["Sample"]
        enzyme = row["Enzyme"].lower()

        match = re.match("^(?P<sample>.+)\.raw$", raw_file_name)
        sample = match.groupdict()["sample"]

        study_meta[sample] = {
            "experiment": experiment,
            "enzyme": enzyme,
            "tissue": tissue_map[(sample_stage,sample_tissue)]
        }

    subsets = list()
    experiments = list()
    enzymes = list()
    tissues = list()
    for _,row in ds_tab.iterrows():
        sample = row["sample"]
        experiment = study_meta[sample]["experiment"]
        enzyme = study_meta[sample]["enzyme"]
        tissue = study_meta[sample]["tissue"]

        # File "Adult_CD4Tcells_Gel_Velos_30_f42.raw" does not
        # yield a Comet output file, so we effectively exclude it.
        subset = pd.NA if sample == "Adult_CD4Tcells_Gel_Velos_30_f42" else enzyme

        subsets.append(subset)
        experiments.append(experiment)
        enzymes.append(enzyme)
        tissues.append(tissue)

    return ds_tab.assign(subset=subsets,experiment=experiments,
                         enzyme=enzymes,tissue=tissues)
