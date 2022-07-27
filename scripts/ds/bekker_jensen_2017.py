import os
import re
import shutil
from tempfile import TemporaryDirectory

import pandas as pd
import requests

from ds.utils import parse_sdrf_key_value_field


def enhance_sample_metadata(ds_tab):

    enzyme_map = {
        "Chymotrypsin": "chymotrypsin",
        "Glutamyl endopeptidase": "glu-c",
        "Lys-C": "lys-c",
        "Trypsin": "trypsin"
    }

    sdrf_col_map = {
        "characteristics[organism part]": "tissue",
        "characteristics[cell line]": "cell_line",
        "comment[cleavage agent details]": "enzyme1",
        "comment [cleavage agent details]": "enzyme2",
        "comment[data file]": "raw_file"
    }

    # These 46 samples are absent from the PXD004452
    # SDRF files, so we store their metadata here.
    absent_meta = dict.fromkeys(
        [
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_1",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_2",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_3",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_4",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_5",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_6",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_7",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_8",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_9",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_10",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_11",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_12",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_13",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_14",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_15",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_16",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_17",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_18",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_19",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_20",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_21",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_22",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_23",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_24",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_25",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_26",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_27",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_28",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_29",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_30",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_31",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_32",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_33",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_34",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_35",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_36",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_37",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_38",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_39",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_40",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_41",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_42",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_43",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_44",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_45_2",
            "20151022_QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac_46"
        ],
        {
            "experiment": "QE5_UPLC10_DBJ_SA_HCT116_REP1_46frac",
            "enzyme": "trypsin",
            "tissue": "colon",
            "cell_line": "HCT116"
        }
    )

    excluded_samples = {
        "20151008_QE5_UPLC10_DBJ_SA_COLON_human_46frac_34_151014141957": "excluded - RAW file conversion error",
        "20151008_QE5_UPLC10_DBJ_SA_LIVER_human_46frac_13_151011001441": "excluded - no spectra searched by Comet",
        "20151008_QE5_UPLC10_DBJ_SA_LIVER_human_46frac_14_151011091039": "excluded - no spectra searched by Comet",
        "20151008_QE5_UPLC10_DBJ_SA_LIVER_human_46frac_18": "excluded - RAW file conversion error"
    }

    base_url = "http://ftp.pride.ebi.ac.uk/pride-archive"
    meta_file_names = [
        "sdrf-celllines.tsv",
        "sdrf-tissues.tsv"
    ]

    study_meta = dict()
    with TemporaryDirectory() as tmp_dir:

        for meta_file_name in meta_file_names:
            meta_file_url = "{}/2017/06/PXD004452/{}".format(base_url,meta_file_name)
            local_meta_file = os.path.join(tmp_dir,meta_file_name)

            with requests.get(meta_file_url,stream=True) as r:
                with open(local_meta_file,"wb") as out_f:
                    shutil.copyfileobj(r.raw,out_f)

            df = pd.read_csv(local_meta_file,sep="\t",
                             usecols=lambda k: k in sdrf_col_map).rename(columns=sdrf_col_map)

            for _,row in df.iterrows():

                # For rows with 2 specified enzymes, we want
                # the latter, so we iterate in reverse order.
                enzyme = None
                for k in ("enzyme2","enzyme1"):
                    d = parse_sdrf_key_value_field(row[k])
                    if d is not None:
                        enzyme = enzyme_map[d["NT"]]
                        break

                if enzyme is None:
                    raise ValueError(
                        "no enzyme specified for sample '{}'".format(sample))

                match = re.match("^(?P<sample>.+)\.raw$",row["raw_file"])
                sample = match["sample"]

                # Ensure samples are unique across the two input tables.
                if sample in study_meta:
                    raise ValueError(
                        "dataset metadata has duplicate sample '{}'".format(sample))

                # Check first with, then without, an additional numerical suffix.
                match = re.match("^[0-9]{8}_(?P<experiment>.+)_[0-9]+_[0-9]+$",sample)
                if not match:
                    match = re.match("^[0-9]{8}_(?P<experiment>.+)_[0-9]+$",sample)
                experiment = match["experiment"]

                study_meta[sample] = {
                    "experiment": experiment,
                    "enzyme": enzyme,
                    "tissue": row["tissue"],
                    "cell_line": row.get("cell_line",default=pd.NA)
                }

    subsets = list()
    experiments = list()
    enzymes = list()
    tissues = list()
    cell_lines = list()
    notes = list()
    for sample in ds_tab["sample"]:

        if sample in study_meta:
            sample_meta = study_meta[sample]
        elif sample in absent_meta:
            sample_meta = absent_meta[sample]
        else:
            raise ValueError(
                "no metadata found for sample '{}'".format(sample))

        if sample in excluded_samples:
            note = excluded_samples[sample]
            subset = pd.NA
        else:
            note = pd.NA
            subset = sample_meta["enzyme"]

        subsets.append(subset)
        experiments.append(sample_meta["experiment"])
        enzymes.append(sample_meta["enzyme"])
        tissues.append(sample_meta["tissue"])
        cell_lines.append(sample_meta["cell_line"])
        notes.append(note)

    ds_tab = ds_tab.assign(subset=subsets,experiment=experiments,enzyme=enzymes,
                           tissue=tissues,cell_line=cell_lines,note=notes)
    return ds_tab.sort_values(by=["experiment","sample"])

