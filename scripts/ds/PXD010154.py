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

    base_url = "http://ftp.pride.ebi.ac.uk/pride/data/archive"
    meta_file_name = "Tissues_Rawfile_list.xlsx"
    meta_file_url = "{}/2019/07/PXD010154/{}".format(base_url,meta_file_name)

    with TemporaryDirectory() as tmp_dir:

        local_meta_file = os.path.join(tmp_dir,meta_file_name)
        with requests.get(meta_file_url,stream=True) as r:
            with open(local_meta_file,"wb") as out_f:
                shutil.copyfileobj(r.raw,out_f)

        df = pd.read_excel(local_meta_file,sheet_name="Sheet1",
                           usecols=["Raw files","Tissue name","Proteases"],
                           engine="openpyxl")

    study_meta = dict()
    synthetic_samples = set()
    for _,row in df.iterrows():
        raw_file_name = row["Raw files"]
        sample_tissue = row["Tissue name"]
        protease = row["Proteases"].lower()

        # In "Tissues_Rawfile_list.xlsx","01697_B01_P018020_S00_N02_R1.raw"
        # is listed erroneously as "01697_B01_P018020_S00_N02_R2.raw".
        if raw_file_name == "01697_B01_P018020_S00_N02_R2.raw":
            raw_file_name = "01697_B01_P018020_S00_N02_R1.raw"

        match = re.match("^(?P<sample>.+)\.raw$", raw_file_name)
        sample = match.groupdict()["sample"]

        if sample_tissue == "-":
            synthetic_samples.add(sample)
            continue

        study_meta[sample] = {
            "enzyme": protease,
            "tissue": tissue_map[sample_tissue]
        }

    ds_tab = ds_tab.loc[~ds_tab["sample"].isin(synthetic_samples)]

    subsets = list()
    enzymes = list()
    tissues = list()
    for _,row in ds_tab.iterrows():
        sample = row["sample"]
        enzyme = study_meta[sample]["enzyme"]
        tissue = study_meta[sample]["tissue"]

        # File "01323_D03_P013562_S00_N20_R1.raw" cannot be converted
        # to mzML format without error, so we effectively exclude it.
        subset = pd.NA if sample == "01323_D03_P013562_S00_N20_R1" else enzyme

        subsets.append(subset)
        enzymes.append(enzyme)
        tissues.append(tissue)

    return ds_tab.assign(subset=subsets,enzyme=enzymes,tissue=tissues)