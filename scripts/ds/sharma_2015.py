import re


def enhance_sample_metadata(ds_tab):

    experiments = list()
    enzymes = list()
    for sample in ds_tab["sample"]:
        enzyme = "trypsin"

        match = re.match("^(?P<experiment>.+)_[0-9]+(?:[A-Z]+)?$",sample)
        if match:
            experiment = match["experiment"]
        else:
            experiment = sample

        experiments.append(experiment)
        enzymes.append(enzyme)

    ds_tab = ds_tab.assign(experiment=experiments,enzyme=enzymes)
    return ds_tab.sort_values(by=["experiment","sample"])
