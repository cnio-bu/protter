import re

import pandas as pd


def enhance_sample_metadata(ds_tab):

    rel_experiments = [
        "PI-401",
        "PI-656",
        "PI-702",
        "PI-714",
        "PreV-012",
        "PreV-013",
        "PreV-030",
        "PreV-064"
    ]

    rel_sample_regex = re.compile(
        "^(?P<experiment>{})_Fr.+_R[12](?:_[0-9]+)?$".format("|".join(rel_experiments))
    )

    excluded_samples = {
        "PI-714_Fr1_R2": "excluded - RAW file conversion error",
        "PI-714_Fr2-3_R2": "excluded - RAW file conversion error",
        "PI-714_Fr4-5_R1": "excluded - RAW file conversion error",
        "PI-714_Fr19-20-21_R2": "excluded - RAW file conversion error",
        "PreV-012_Fr14-15_R1": "excluded - RAW file conversion error",
        "PreV-012_FrE-F_R2": "excluded - RAW file conversion error",
        "PreV-030_Fr16-17-18_R2": "excluded - RAW file conversion error"
    }

    rel_mask = list()
    experiments = list()
    for sample in ds_tab["sample"]:
        match = rel_sample_regex.match(sample)
        if match:
            experiments.append(match["experiment"])
            is_rel_sample = True
        else:
            is_rel_sample = False
        rel_mask.append(is_rel_sample)

    ds_tab = ds_tab.loc[rel_mask]

    subsets = list()
    enzymes = list()
    notes = list()
    for sample in ds_tab["sample"]:
        enzyme = "trypsin"

        if sample in excluded_samples:
            note = excluded_samples[sample]
            subset = pd.NA
        else:
            note = pd.NA
            subset = enzyme

        subsets.append(subset)
        enzymes.append(enzyme)
        notes.append(note)

    ds_tab = ds_tab.assign(subset=subsets,experiment=experiments,
                           enzyme=enzymes,note=notes)
    return ds_tab.sort_values(by=["experiment","sample"])
