#!/usr/bin/env python3
"""Prepare workflow sample sheet file.

This is primarily intended to create an initial sample
sheet file, but you should always review the file before
using it, and you may need to make further changes.
"""

from argparse import ArgumentParser

import pandas as pd

from common import get_dataset_metadata,load_config_file
from ds import dataset_module_names,load_dataset_module


def init_sample_metadata(ds,config):
    ds_meta = get_dataset_metadata(ds,config)

    recs = list()
    for sample in sorted(ds_meta.keys()):
        sample_meta = ds_meta[sample]
        recs.append({
            "dataset": ds,
            "project": sample_meta.get("project",pd.NA),
            "sample": sample,
            "file": sample_meta["file"],
            "checksum": sample_meta.get("checksum",pd.NA)
        })

    col_names = ["dataset","project","sample","file","checksum"]
    ds_tab = pd.DataFrame(recs,columns=col_names,dtype="string")
    ds_tab.dropna(axis=1,how="all",inplace=True)

    return ds_tab


if __name__ == "__main__":

    ap = ArgumentParser(description=__doc__)
    ap.add_argument("config_file",
                    help="input config YAML file")
    ap.add_argument("sample_file",
                    help="output sample sheet TSV file")
    args = ap.parse_args()

    config_file = args.config_file
    sample_file = args.sample_file


    config = load_config_file(config_file)

    datasets = [ds for ds in config["datasets"]
                if config["datasets"][ds].get("enabled",False)]
    known_datasets = dataset_module_names()

    ds_tabs = list()
    subsetted_datasets = set()
    for ds in datasets:
        ds_tab = init_sample_metadata(ds,config)

        if ds in known_datasets and len(ds_tab) > 0:
            ds_module = load_dataset_module(ds)
            ds_tab = ds_module.enhance_sample_metadata(ds_tab)

        if "subset" in ds_tab.columns:
            subsetted_datasets.add(ds)

        ds_tabs.append(ds_tab)

    if subsetted_datasets:
        for i,ds in enumerate(datasets):
            if ds not in subsetted_datasets:
                ds_tabs[i] = ds_tabs[i].assign(subset=["all"] * len(ds_tabs[i]))

    samples = pd.concat(ds_tabs)
    samples.to_csv(sample_file,sep="\t",na_rep="NA",
                   index=False,encoding="utf-8")
