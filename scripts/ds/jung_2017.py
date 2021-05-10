import re

import pandas as pd


def enhance_sample_metadata(ds_tab):
    ds_tab["project"] = ["MSV000079789"] * len(ds_tab)

    excluded_samples = {
        "M0209PA_9758_1_V_21": "excluded - RAW file conversion error"
    }

    subsets = list()
    experiments = list()
    enzymes = list()
    notes = list()
    for sample in ds_tab["sample"]:
        enzyme = "trypsin"

        match = re.match("^(?P<experiment>.+)_[0-9]+(?:\+[0-9]+)?$",sample)
        experiment = match["experiment"]

        if sample in excluded_samples:
            note = excluded_samples[sample]
            subset = pd.NA
        else:
            note = pd.NA
            subset = enzyme

        subsets.append(subset)
        experiments.append(experiment)
        enzymes.append(enzyme)
        notes.append(note)

    ds_tab = ds_tab.assign(subset=subsets,experiment=experiments,
                           enzyme=enzymes,note=notes)
    return ds_tab.sort_values(by=["experiment","sample"])