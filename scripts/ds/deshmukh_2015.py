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
        "Characteristics[cell line]": "cell_line",
        "comment[cleavage agent details]": "enzyme",
        "comment[data file]": "raw_file",
        "sampleAccession": "biosample"
    }

    tissue_map = {
        "Triceps muscle": "triceps"
    }

    base_url = "http://ftp.pride.ebi.ac.uk/pride/data/archive"
    meta_file_name = "sdrf.tsv"
    meta_file_url = "{}/2015/04/PXD000288/{}".format(base_url,meta_file_name)

    with TemporaryDirectory() as tmp_dir:
        local_meta_file = os.path.join(tmp_dir,meta_file_name)
        with requests.get(meta_file_url,stream=True) as r:
            with open(local_meta_file,"wb") as out_f:
                shutil.copyfileobj(r.raw,out_f)

        df = pd.read_csv(local_meta_file,sep="\t",usecols=sdrf_col_map.keys(),
                         na_values=["not applicable","not available"],
                         keep_default_na=False).convert_dtypes().rename(columns=sdrf_col_map)

    study_meta = dict()
    for _,row in df.iterrows():

        tissue = row["tissue"]
        if not pd.isna(tissue):
            tissue = tissue_map[tissue]

        d = parse_sdrf_key_value_field(row["enzyme"])
        enzyme = enzyme_map[d["NT"]]

        match = re.match("^(?P<sample>.+)\.raw$",row["raw_file"])
        sample = match["sample"]

        study_meta[sample] = {
            "biosample": row["biosample"],
            "enzyme": enzyme,
            "tissue": tissue,
            "cell_line": row["cell_line"]
        }

    biosamples = list()
    enzymes = list()
    tissues = list()
    cell_lines = list()
    for sample in ds_tab["sample"]:
        biosamples.append(study_meta[sample]["biosample"])
        enzymes.append(study_meta[sample]["enzyme"])
        tissues.append(study_meta[sample]["tissue"])
        cell_lines.append(study_meta[sample]["cell_line"])

    ds_tab = ds_tab.assign(biosample=biosamples,enzyme=enzymes,
                           tissue=tissues,cell_line=cell_lines)
    return ds_tab.sort_values(by=["biosample","sample"])
