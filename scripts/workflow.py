from calendar import timegm
from collections import defaultdict
import json
import os
import shutil
import time

import pandas as pd
import snakemake

from .common import (dataset_dir,
                     dataset_source,
                     get_dataset_metadata,
                     get_pride_dataset_readme_url,
                     get_pride_file_mtime_info,
                     get_samples,
                     is_comet_fmt,
                     is_msconvert_fmt,
                     url_basename,
                     utc_strptime)


def comet_input_file(wildcards,config,samples):
    ds = wildcards.ds
    ds_fmt = config["datasets"][ds]["fmt"]

    if is_comet_fmt(ds_fmt,config):
        input_file = sample_data_file(wildcards,config,samples)
    elif is_msconvert_fmt(ds_fmt,config):
        input_file = msconvert_output_pattern(config)
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


def dataset_metadata(ds,subset,config,samples):
    ds_dir = dataset_dir(ds,config)
    ds_meta_file = os.path.join(ds_dir,"dataset-metadata.json")
    with open(ds_meta_file) as f:
        ds_meta = json.load(f)
    ds_samples = samples.xs(key=ds,level="dataset")
    if subset == "all":
        subset_sample_ids = list(ds_samples["sample"])
    elif "subset" in ds_samples.columns:
        subset_sample_ids = list(ds_samples[ds_samples["subset"] == subset,"samples"])
    else:
        raise ValueError("subset not configured: '{}'".format(subset))
    subset_meta = {}
    for sample_id in subset_sample_ids:
        try:
            subset_meta[sample_id] = ds_meta["samples"][sample_id]
        except KeyError:
            raise ValueError("unknown sample: '{}'".format(sample_id))
    ds_meta["samples"] = subset_meta
    return ds_meta


def dataset_subsets(ds,samples):
    ds_samples = samples.xs(key=ds,level="dataset")
    if "subset" in ds_samples.columns:
        subsets = list(ds_samples["subset"].dropna().unique())
        if not subsets:
            subsets = ["all"]
    else:
        subsets = ["all"]
    return subsets


def download_sample_output_pattern(config):
    out_patt = os.path.join(config["dataset_path"],"{ds}","{sample}.{fmt}{gzip_ext}")

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


def msconvert_input_file(wildcards,config,samples):
    ds = wildcards.ds
    ds_fmt = config["datasets"][ds]["fmt"]

    if is_msconvert_fmt(ds_fmt,config):
        msconvert_input = sample_data_file(wildcards,config,samples)
    elif is_comet_fmt(ds_fmt,config):  # i.e. no need to run msconvert
        msconvert_input = ""
    else:
        raise ValueError(
            "unsupported proteomics file format: '{}'".format(ds_fmt))

    return msconvert_input


def msconvert_output_pattern(config):
    out_patt = os.path.join(config["dataset_path"],"{ds}","{sample}.mzML.gz")

    try:
        output_flags = config["rules"]["msconvert"]["output"]["flags"]
    except KeyError:
        output_flags = []

    # Call the flag functions right-to-left, as
    # if in a series of nested function calls.
    for flag_name in reversed(output_flags):
        flag_func = getattr(snakemake.io, flag_name)
        out_patt = flag_func(out_patt)

    return out_patt


def msconvert_rule_config(config):
    raw_config = config["software"]["msconvert"]
    rule_config = defaultdict(str)

    rule_config["docker_image"] = raw_config["docker_image"]
    rule_config["docker_host"] = raw_config.get("docker_host","local")
    if rule_config["docker_host"] != "local":
        remote_host = raw_config["remote"]["hostname"]
        remote_user = raw_config["remote"]["username"]
        rule_config["remote_addr"] = "{}@{}".format(remote_user,remote_host)
        rule_config["remote_key"] = raw_config["remote"].get("private_key","")

        local_host = raw_config["local"]["hostname"]
        local_user = raw_config["local"]["username"]
        rule_config["local_addr"] = "{}@{}".format(local_user,local_host)
        rule_config["local_key"] = raw_config["local"].get("private_key","")

    return rule_config


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
    ds_meta = dataset_metadata(ds,"all",config,samples)
    sample_meta = ds_meta["samples"][sample]
    if "url" in sample_meta:
        file_name = url_basename(sample_meta["url"])
        ds_dir = dataset_dir(ds,config)
        input_file = os.path.join(ds_dir,file_name)
    elif "path" in sample_meta:
        input_file = sample_meta["path"]
    else:
        raise ValueError(
            "cannot get local file for sample: '{}'".format(sample))
    return input_file


def sync_dataset_metadata(ds,config):
    ds_dir = dataset_dir(ds,config)
    ds_meta_file = os.path.join(ds_dir,"dataset-metadata.json")
    try:
        with open(ds_meta_file) as f:
            ds_meta = json.load(f)
    except FileNotFoundError:
        update_needed = True
    else:
        if ds_meta["config"] != config["datasets"][ds]:
            update_needed = True
        else:
            meta_file_mtime = os.path.getmtime(ds_meta_file)
            ds_src = dataset_source(ds,config)
            if ds_src == "local":
                latest_mtime = max([os.path.getmtime(os.path.join(ds_dir,x))
                                    for x in os.listdir(ds_dir)])
            elif ds_src == "PRIDE":
                # Take README 'last-modified' time as representative,
                # as it contains metadata on the other dataset files.
                readme_url = get_pride_dataset_readme_url(ds,ds_meta)
                mtime_info = get_pride_file_mtime_info([readme_url])
                latest_mtime = mtime_info[readme_url]
            else:
                raise ValueError(
                    "unsupported proteomics data source: '{}'".format(ds_src))

            update_needed = meta_file_mtime < latest_mtime

    if update_needed:
        tmp_ds_meta_file = os.path.join(ds_dir,"dataset-metadata-tmp.json")
        os.makedirs(ds_dir,exist_ok=True)
        try:
            with open(tmp_ds_meta_file,"x") as f:
                ds_meta = get_dataset_metadata(ds,config)
                json.dump(ds_meta,f)
        except FileExistsError:
            timeout = 600
            stop_time = time.time() + timeout
            while not os.path.isfile(ds_meta_file):
                if time.time() >= stop_time:
                    raise TimeoutError(
                        "waited {} seconds for creation of metadata"
                        " file for dataset '{}'".format(timeout,ds))
            with open(ds_meta_file) as f:
                ds_meta = json.load(f)
            if ds_meta["config"] != config["datasets"][ds]:
                raise ValueError("")
        else:
            shutil.move(tmp_ds_meta_file,ds_meta_file)


def sync_sample_proxy_files(ds,config,samples):

    ds_dir = dataset_dir(ds,config)
    marked_proxy_files = set()
    if os.path.exists(ds_dir):
        for item_name in os.listdir(ds_dir):
            item_path = os.path.join(ds_dir,item_name)
            if os.path.isfile(item_path) and item_path.endswith("_proxy.json"):
                marked_proxy_files.add(item_path)

    if dataset_source(ds,config) != "local":
        os.makedirs(ds_dir,exist_ok=True)
        ds_meta = dataset_metadata(ds,"all",config,samples)
        for sample,sample_meta in ds_meta["samples"].items():
            file_url = sample_meta["url"]
            mod_dt = utc_strptime(sample_meta["last-modified"])
            file_mtime = timegm(mod_dt.timetuple())
            proxy_file_name = "{}_proxy.json".format(url_basename(file_url))
            proxy_file_path = os.path.join(ds_dir,proxy_file_name)
            with open(proxy_file_path,"w") as f:
                json.dump(sample_meta,f)
            os.utime(proxy_file_path,(time.time(),file_mtime))
            marked_proxy_files.discard(proxy_file_path)

    for marked_proxy_file in marked_proxy_files:
        os.remove(marked_proxy_file)
