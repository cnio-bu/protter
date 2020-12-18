#!/usr/bin/env python3
"""Initialise workflow sample sheet file."""

from argparse import ArgumentParser
import csv

from common import get_dataset_metadata,load_config_file


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

    sample_config_recs = list()
    for ds in config["datasets"]:
        if not config["datasets"][ds].get("enabled",False):
            continue

        ds_meta = get_dataset_metadata(ds,config)

        for sample in sorted(ds_meta["samples"].keys()):
            sample_config_recs.append({
                "dataset": ds,
                "sample": sample
            })

    col_names = ["dataset","sample"]
    with open(sample_file,"w",newline="") as f:
        writer = csv.DictWriter(f,col_names,delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for sample_config_rec in sample_config_recs:
            writer.writerow(sample_config_rec)
