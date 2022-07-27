import os
import re
import shutil
from tempfile import TemporaryDirectory

import pandas as pd
import requests


def enhance_sample_metadata(ds_tab):

    tissue_map = {
        "adrenal gland": "adrenal gland",
        "appendix/vermiform appendix": "appendix",
        "bone marrow": "bone marrow",
        "brain": "brain",
        "colon": "colon",
        "duodenum": "duodenum",
        "endometrium/uterine endometrium": "endometrium",
        "esophagus": "esophagus",
        "fallopian tube/oviduct": "fallopian tube",
        "fat/adipose tissue": "fat",
        "gallbladder": "gallbladder",
        "heart": "heart",
        "kidney": "kidney",
        "liver": "liver",
        "lung": "lung",
        "lymph node": "lymph node",
        "ovary": "ovary",
        "pancreas": "pancreas",
        "pituitary/hypophysis": "pituitary gland",
        "placenta": "placenta",
        "prostate": "prostate",
        "rectum": "rectum",
        "salivary gland": "salivary gland",
        "small intestine": "small intestine",
        "smooth muscle": "smooth muscle",
        "spleen": "spleen",
        "stomach": "stomach",
        "testis": "testis",
        "thyroid": "thyroid",
        "tonsil": "tonsil",
        "urinary bladder": "urinary bladder"
    }

    excluded_samples = {
        "01323_D03_P013562_S00_N20_R1": "excluded - RAW file conversion error"
    }

    base_url = "http://ftp.pride.ebi.ac.uk/pride-archive"
    meta_file_name = "Tissues_Rawfile_list.xlsx"
    meta_file_url = "{}/2019/07/PXD010154/{}".format(base_url,meta_file_name)

    with TemporaryDirectory() as tmp_dir:

        local_meta_file = os.path.join(tmp_dir,meta_file_name)
        with requests.get(meta_file_url,stream=True) as r:
            with open(local_meta_file,"wb") as out_f:
                shutil.copyfileobj(r.raw,out_f)

        df = pd.read_excel(local_meta_file,sheet_name="Sheet1",
                           usecols=["Raw files","Sample source","Tissue name",
                                    "Experiments","Proteases"],
                           engine="openpyxl")

    study_meta = dict()
    synthetic_samples = set()
    for _,row in df.iterrows():
        raw_file_name = row["Raw files"].strip()
        sample_source = row["Sample source"].strip()
        experiment = row["Experiments"].strip()
        sample_tissue = row["Tissue name"].strip()
        protease = row["Proteases"].lower()

        # In "Tissues_Rawfile_list.xlsx","01697_B01_P018020_S00_N02_R1.raw"
        # is listed erroneously as "01697_B01_P018020_S00_N02_R2.raw".
        if raw_file_name == "01697_B01_P018020_S00_N02_R2.raw":
            raw_file_name = "01697_B01_P018020_S00_N02_R1.raw"

        match = re.match("^(?P<sample>.+)\.raw$", raw_file_name)
        sample = match["sample"]

        if sample_source == "Synthetic":
            synthetic_samples.add(sample)
            continue

        study_meta[sample] = {
            "experiment": experiment,
            "enzyme": protease,
            "tissue": tissue_map[sample_tissue]
        }

    ds_tab = ds_tab.loc[~ds_tab["sample"].isin(synthetic_samples)]

    subsets = list()
    experiments = list()
    enzymes = list()
    tissues = list()
    notes = list()
    for sample in ds_tab["sample"]:
        experiment = study_meta[sample]["experiment"]
        enzyme = study_meta[sample]["enzyme"]
        tissue = study_meta[sample]["tissue"]

        if sample in excluded_samples:
            note = excluded_samples[sample]
            subset = pd.NA
        else:
            note = pd.NA
            subset = enzyme

        subsets.append(subset)
        experiments.append(experiment)
        enzymes.append(enzyme)
        tissues.append(tissue)
        notes.append(note)

    ds_tab = ds_tab.assign(subset=subsets,experiment=experiments,
                           enzyme=enzymes,tissue=tissues,notes=notes)
    return ds_tab.sort_values(by=["experiment","sample"])
