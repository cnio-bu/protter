import re


def enhance_sample_metadata(ds_tab):
    enzymes = ["trypsin"] * len(ds_tab)

    experiments = list()
    for _,row in ds_tab.iterrows():
        match = re.match("^(?P<experiment>.+)_Fr[0-9]+$",row["sample"])
        experiment = match.groupdict()["experiment"]
        experiments.append(experiment)

    return ds_tab.assign(experiment=experiments,enzyme=enzymes)
