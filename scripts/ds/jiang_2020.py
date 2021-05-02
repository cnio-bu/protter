import re


def enhance_sample_metadata(ds_tab):
    enzymes = ["trypsin"] * len(ds_tab)

    experiments = list()
    for _,row in ds_tab.iterrows():
        match = re.match("^(?P<experiment>.+)_Fr[0-9]+$",row["sample"])
        experiments.append(match["experiment"])

    ds_tab = ds_tab.assign(experiment=experiments,enzyme=enzymes)
    return ds_tab.sort_values(by=["experiment","sample"])
