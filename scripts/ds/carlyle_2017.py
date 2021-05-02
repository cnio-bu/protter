import os
import re
import shutil
from tempfile import TemporaryDirectory

import pandas as pd
import requests

from ds.utils import parse_sdrf_key_value_field


def enhance_sample_metadata(ds_tab):

    enzyme_map = {
        "Trypsin": "trypsin"
    }

    sdrf_col_map = {
        "characteristics[organism part]": "tissue",
        "comment[cleavage agent details]": "enzyme",
        "comment[data file]": "raw_file",
        "sampleAccession": "biosample"
    }

    base_url = "http://ftp.pride.ebi.ac.uk/pride/data/archive"
    meta_file_name = "sdrf.tsv"
    meta_file_url = "{}/2018/01/PXD005445/{}".format(base_url,meta_file_name)

    with TemporaryDirectory() as tmp_dir:

        local_meta_file = os.path.join(tmp_dir,meta_file_name)
        with requests.get(meta_file_url,stream=True) as r:
            with open(local_meta_file,"wb") as out_f:
                shutil.copyfileobj(r.raw,out_f)

        df = pd.read_csv(local_meta_file,sep="\t",
                         usecols=sdrf_col_map.keys()).rename(columns=sdrf_col_map)

    study_meta = dict()
    control_samples = set()
    for _,row in df.iterrows():
        biosample = row["biosample"]
        tissue = row["tissue"]

        d = parse_sdrf_key_value_field(row["enzyme"])
        enzyme = enzyme_map[d["NT"]]

        match = re.match("^(?P<sample>.+)\.raw$",row["raw_file"])
        sample = match["sample"]

        if tissue in ("not applicable","not available"):
            control_samples.add(sample)
            continue

        study_meta[sample] = {
            "biosample": biosample,
            "enzyme": enzyme,
            "tissue": tissue
        }

    ds_tab = ds_tab.loc[~ds_tab["sample"].isin(control_samples)]

    biosamples = list()
    enzymes = list()
    tissues = list()
    for _,row in ds_tab.iterrows():
        sample = row["sample"]

        biosamples.append(study_meta[sample]["biosample"])
        enzymes.append(study_meta[sample]["enzyme"])
        tissues.append(study_meta[sample]["tissue"])

    ds_tab = ds_tab.assign(biosample=biosamples,enzyme=enzymes,tissue=tissues)
    return ds_tab.sort_values(by=["biosample","sample"])
