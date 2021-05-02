import re


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

    experiments = list()
    dropped_samples = set()
    for _,row in ds_tab.iterrows():
        sample = row["sample"]

        match = rel_sample_regex.match(sample)
        if match:
            experiments.append(match["experiment"])
        else:
            dropped_samples.add(sample)

    ds_tab = ds_tab.loc[~ds_tab["sample"].isin(dropped_samples)]
    enzymes = ["trypsin"] * len(ds_tab)

    ds_tab = ds_tab.assign(experiment=experiments,enzyme=enzymes)
    return ds_tab.sort_values(by=["experiment","sample"])
