import re


def enhance_sample_metadata(ds_tab):
    enzymes = ["trypsin"] * len(ds_tab)

    groups = list()
    for _,row in ds_tab.iterrows():
        match = re.match("^(?P<group>.+)_Fr[0-9]+$",row["sample"])
        group = match.groupdict()["group"]
        groups.append(group)

    return ds_tab.assign(group=groups,enzyme=enzymes)
