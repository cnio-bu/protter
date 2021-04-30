import os
import re
import shutil
from tempfile import TemporaryDirectory

import pandas as pd
import requests


def enhance_sample_metadata(ds_tab):
    enzymes = ["trypsin"] * len(ds_tab)

    tissue_map = {
        "1 - right colon": "right colon",
        "2 - transverse colon": "transverse colon",
        "3 - left colon": "left colon",
        "4 - rectum": "rectum"
    }

    # We use supplementary information file for Zhang et al. (2014), because
    # not all of this study's PRIDE projects have associated 'sdrf.tsv' files.
    base_url = "https://static-content.springer.com/esm/art%3A10.1038%2Fnature13438/MediaObjects"
    meta_file_name = "41586_2014_BFnature13438_MOESM5_ESM.xlsx"
    meta_file_url = "{}/{}".format(base_url,meta_file_name)

    with TemporaryDirectory() as tmp_dir:
        local_meta_file = os.path.join(tmp_dir,meta_file_name)
        with requests.get(meta_file_url,stream=True) as r:
            with open(local_meta_file,"wb") as out_f:
                shutil.copyfileobj(r.raw,out_f)

        study_meta = pd.read_excel(local_meta_file,sheet_name="S1_sample_annotation",
                                   usecols=["CPTAC sample ID","Tumor.site"],skiprows=56,
                                   engine="openpyxl")
        study_meta.set_index("CPTAC sample ID",drop=True,inplace=True,verify_integrity=True)
        tumor_sites = study_meta["Tumor.site"]

    biosample_regex = re.compile("^(?P<biosample>TCGA-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+)",
                                 re.IGNORECASE)

    tissues = list()
    biosamples = list()
    for _,row in ds_tab.iterrows():
        match = biosample_regex.match(row["sample"])
        biosample = match["biosample"]

        tumor_site = tumor_sites[biosample]
        tissue = tissue_map[tumor_site]

        biosamples.append(biosample)
        tissues.append(tissue)

    ds_tab = ds_tab.assign(biosample=biosamples,enzyme=enzymes,tissue=tissues)
    return ds_tab.sort_values(by=["biosample","sample"])
