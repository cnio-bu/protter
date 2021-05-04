import re


def enhance_sample_metadata(ds_tab):
    enzymes = ["trypsin"] * len(ds_tab)

    experiments = list()
    for sample in ds_tab["sample"]:
        match = re.match("^(?P<experiment>.+)_Fr[0-9]+$",sample)
        experiments.append(match["experiment"])

    ds_tab = ds_tab.assign(experiment=experiments,enzyme=enzymes)
    return ds_tab.sort_values(by=["experiment","sample"])
