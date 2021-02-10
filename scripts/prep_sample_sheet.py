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
    recs = list()
    ds_meta = get_dataset_metadata(ds,config)
    for sample in sorted(ds_meta["samples"].keys()):
        recs.append({
            "dataset": ds,
            "sample": sample
        })
    col_names = ["dataset","sample"]
    return pd.DataFrame(recs,columns=col_names,dtype="string")


if __name__ == "__main__":

    ap = ArgumentParser(description=__doc__)
    ap.add_argument("--config-file",
                    help="input config YAML file")
    ap.add_argument("--sample-file",
                    help="output sample sheet TSV file")
    args = ap.parse_args()

    config_file = args.config_file
    sample_file = args.sample_file


    config = load_config_file(config_file)

    datasets = [ds for ds in config["datasets"]
                if config["datasets"][ds].get("enabled",False)]
    known_datasets = dataset_module_names()

    ds_tabs = list()
    for ds in datasets:
        ds_tab = init_sample_metadata(ds,config)
        if ds in known_datasets:
            ds_module = load_dataset_module(ds)
            ds_tab = ds_module.enhance_sample_metadata(ds_tab)
        ds_tabs.append(ds_tab)

    samples = pd.concat(ds_tabs)
    samples.to_csv(sample_file,sep="\t",na_rep="NA",
                   index=False,encoding="utf-8")
