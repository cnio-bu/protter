from collections import defaultdict
import os

import pandas as pd
import snakemake

from .common import (dataset_dir,
                     get_samples,
                     is_comet_fmt,
                     is_raw_fmt,
                     is_wget_url,
                     url_basename)


def comet_input_file(wildcards,config,samples):
    ds = wildcards.ds
    ds_fmt = config["datasets"][ds]["fmt"]

    if is_comet_fmt(ds_fmt,config):
        input_file = sample_data_file(wildcards,config,samples)
    elif is_raw_fmt(ds_fmt,config):
        input_file = mzml_output_pattern(config)
    else:
        raise ValueError(
            "unsupported proteomics file format: '{}'".format(ds_fmt))
    return input_file


def dataset_dbs(ds,config):
    return [db for db in config["dbs"]
            if config["dbs"][db]["enabled"]
            and db in config["datasets"][ds]["dbs"]]


def dataset_groupings(ds,config):
    ds_conf = config["datasets"][ds]
    if "groupings" in ds_conf:
        groupings = ds_conf["groupings"]
    else:
        groupings = config["default_groupings"]
    return groupings


def dataset_subsets(ds,samples):
    ds_samples = samples.xs(key=ds,level="dataset")
    if "subset" in ds_samples.columns:
        subsets = list(ds_samples["subset"].dropna().unique())
        if not subsets:
            subsets = ["all"]
    else:
        subsets = ["all"]
    return subsets


def download_sample_file_checksum(wildcards,downloads):
    checksum = downloads.loc[(wildcards.ds,wildcards.dl_file),"checksum"]
    if pd.isna(checksum):
        checksum = ""
    return checksum


def download_sample_file_url(wildcards,downloads):
    return downloads.loc[(wildcards.ds,wildcards.dl_file),"file"]


def download_sample_output_pattern(config):
    out_patt = os.path.join(config["dataset_path"],"{ds}","{dl_file}")

    try:
        output_flags = config["rules"]["download_sample"]["output"]["flags"]
    except KeyError:
        output_flags = []

    # Call the flag functions right-to-left, as
    # if in a series of nested function calls.
    for flag_name in reversed(output_flags):
        flag_func = getattr(snakemake.io, flag_name)
        out_patt = flag_func(out_patt)

    return out_patt


def raw_input_file(wildcards,config,samples):
    ds = wildcards.ds
    ds_fmt = config["datasets"][ds]["fmt"]

    if is_raw_fmt(ds_fmt,config):
        raw_input = sample_data_file(wildcards,config,samples)
    elif is_comet_fmt(ds_fmt,config):  # i.e. no need to convert
        raw_input = ""
    else:
        raise ValueError(
            "unsupported proteomics file format: '{}'".format(ds_fmt))

    return raw_input


def mzml_output_pattern(config):
    out_patt = os.path.join(config["dataset_path"],"{ds}","{sample}.mzML.gz")

    try:
        output_flags = config["rules"]["raw_to_mzml"]["output"]["flags"]
    except KeyError:
        output_flags = []

    # Call the flag functions right-to-left, as
    # if in a series of nested function calls.
    for flag_name in reversed(output_flags):
        flag_func = getattr(snakemake.io, flag_name)
        out_patt = flag_func(out_patt)

    return out_patt


def percolator_enzyme(wildcards,config,samples):
    ds = wildcards.ds
    subset = wildcards.subset
    grouping = wildcards.grouping
    group = wildcards.group

    rel_samples = samples.xs(key=ds,level="dataset")

    if subset != "all":
        if "subset" not in rel_samples.columns:
            raise ValueError("subset not configured: '{}'".format(subset))
        rel_samples = rel_samples.loc[rel_samples["subset"] == subset]

    if "enzyme" in rel_samples.columns:
        if grouping in rel_samples.columns:
            group_enzymes = set(rel_samples.loc[rel_samples[grouping] == group,
                                                "enzyme"].dropna())
            if len(group_enzymes) != 1:
                raise ValueError(
                    "enzyme not correctly configured for group: '{}'".format(group))
            enzyme = group_enzymes.pop()
        elif grouping == "single":  # i.e. group == sample
            enzyme = rel_samples.loc[group,"enzyme"]
            if pd.isna(enzyme):
                raise ValueError(
                    "enzyme not correctly configured for sample: '{}'".format(group))
        else:
            raise ValueError("grouping not configured: '{}'".format(grouping))
    else:
        enzyme = config["software"]["percolator"]["default_enzyme"]

    if enzyme not in config["software"]["percolator"]["enzymes"]:
        raise ValueError("unknown enzyme: '{}'".format(enzyme))
    return enzyme


def percolator_input_files(wc,samples):
    group_to_samples = get_samples(wc.ds,wc.subset,wc.grouping,samples)
    return ["out/{db}/split_pins/{ds}/{subset}/{sample}.{sdb}.pin".format(
            db=wc.db,ds=wc.ds,subset=wc.subset,sample=sample,sdb=wc.sdb)
            for sample in group_to_samples[wc.group]]


def sample_data_file(wildcards,config,samples):
    ds = wildcards.ds
    sample = wildcards.sample
    input_file = samples.loc[(ds,sample),"file"]

    if is_wget_url(input_file):
        file_name = url_basename(input_file)
        ds_dir = dataset_dir(ds,config)
        input_file = os.path.join(ds_dir,file_name)

    return input_file
