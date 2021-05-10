import re

import pandas as pd


def enhance_sample_metadata(ds_tab):

    rel_experiments = [
        "Swovick_M2D",
        "Swovick_2M2D",
        "Swovick_M4D",
        "Swovick_2M4D",
        "Swovick_M6D",
        "Swovick_2M6D"
    ]
    rel_exp_patt = "|".join(re.escape(x) for x in rel_experiments)

    rel_sample_regex = re.compile(
        "^(?P<experiment>{})-E?[1-8]_[0-9]+(?:-[0-9]+)+$".format(rel_exp_patt)
    )

    excluded_samples = {
        "Swovick_M4D-1_17-037": "excluded - RAW file conversion error"
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
