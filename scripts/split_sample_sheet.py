#!/usr/bin/env python3
"""Split sample sheet file."""

from argparse import ArgumentDefaultsHelpFormatter,ArgumentParser

import pandas as pd


if __name__ == "__main__":

    ap = ArgumentParser(description=__doc__,
                        formatter_class=ArgumentDefaultsHelpFormatter)
    ap.add_argument("input_file",
                    help="Input sample sheet TSV file (e.g. 'samples.tsv').")
    ap.add_argument("output_prefix",
                    help="Output prefix of split sample sheet TSV files (e.g. 'samples_').")
    ap.add_argument("--split-by",choices=["dataset","subset"],default="dataset",
                    help="Key to use when splitting sample sheet. Splitting by dataset makes"
                         " an output sample sheet for each dataset. Splitting by subset makes"
                         " an output sample sheet for each subset of each dataset.")
    args = ap.parse_args()

    input_file = args.input_file
    output_prefix = args.output_prefix
    split_by = args.split_by

    samples = pd.read_csv(input_file,sep="\t",na_values="NA",
                          keep_default_na=False,encoding="utf-8")

    # Whether splitting by dataset or subset,
    # we always need to split by dataset.
    split_samples = {k: x for k,x in samples.groupby("dataset")}

    if split_by == "subset":
        samples_by_subset = dict()
        for ds,ds_samples in split_samples.items():

            ds_samples = ds_samples.loc[ds_samples["subset"].notna()]

            for subset,subset_samples in ds_samples.groupby("subset"):
                split_key = "{}_{}".format(ds,subset)
                samples_by_subset[split_key] = subset_samples

        split_samples = samples_by_subset

    for split_key,sample_sheet in split_samples.items():
        output_file = "{}{}.tsv".format(output_prefix,split_key)
        sample_sheet.dropna(axis=1,how="all",inplace=True)
        sample_sheet.to_csv(output_file,sep="\t",na_rep="NA",
                            index=False,encoding="utf-8")
