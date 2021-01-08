from calendar import timegm
from collections import defaultdict
import json
import os
import time

import snakemake

from .common import (dataset_dir,
                     dataset_source,
                     get_dataset_metadata,
                     get_pride_dataset_readme_url,
                     get_pride_file_mtime_info,
                     is_comet_fmt,
                     is_msconvert_fmt,
                     split_gzip_ext,
                     url_basename,
                     utc_strptime)


def _iter_list_file(file_path):
    with open(file_path) as f:
        for line in f:
            contents = line.strip()
            if contents:
                yield contents


def comet_input_file(wildcards,config):
    ds = wildcards.ds
    ds_fmt = config["datasets"][ds]["fmt"]

    if is_comet_fmt(ds_fmt,config):
        input_file = sample_data_file(wildcards,config)
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
        groupings = config["grouping_default"]
    return groupings


def dataset_metadata(ds,subset,config):
    ds_dir = dataset_dir(ds,config)
    ds_meta_file = os.path.join(ds_dir,"dataset-metadata.json")
    with open(ds_meta_file) as f:
        ds_meta = json.load(f)
    if subset != "all":
        try:
            ds_conf = config["datasets"][ds]
        except KeyError:
            raise ValueError("dataset not configured: '{}'".format(ds))
        try:
            subset_conf = ds_conf["subsets"][subset]
        except KeyError:
            raise ValueError("subset not configured: '{}'".format(subset))
        try:
            sample_file = subset_conf["samples"]
        except KeyError:
            raise ValueError(
                "sample file not configured for subset: '{}'".format(subset))
        subset_meta = {}
        for sample in _iter_list_file(sample_file):
            try:
                subset_meta[sample] = ds_meta["samples"][sample]
            except KeyError:
                raise ValueError("unknown sample: '{}'".format(sample))
        ds_meta["samples"] = subset_meta
    return ds_meta


def dataset_subsets(ds,config):
    if "subsets" in config["datasets"][ds]:
        subsets = list(config["datasets"][ds]["subsets"].keys())
        if "all" in subsets and "samples" in subsets["all"]:
            raise ValueError("cannot configure samples for subset 'all'")
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


def get_samples(ds,subset,grouping,config):
    #get the grouping statement from the config file and create a "gf" function with it
    exec('gf = lambda x: {}'.format(dataset_groupings(ds,config)[grouping]), globals())
    groups = {}
    ds_meta = dataset_metadata(ds,subset,config)
    for f in ds_meta["samples"].keys():
        g = gf(f)
        try:
            groups[g].append(f)
        except KeyError:
            groups[g] = [f]
    return groups


def msconvert_input_file(wildcards,config):
    ds = wildcards.ds
    ds_fmt = config["datasets"][ds]["fmt"]

    if is_msconvert_fmt(ds_fmt,config):
        msconvert_input = sample_data_file(wildcards,config)
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


def percolator_enzyme(wildcards,config):
    ds = wildcards.ds
    subset = wildcards.subset
    try:
        ds_conf = config["datasets"][ds]
    except KeyError:
        raise ValueError("dataset not configured: '{}'".format(ds))
    default_enzyme = config["software"]["percolator"]["default_enzyme"]
    try:
        subset_conf = ds_conf["subsets"][subset]
    except KeyError:
        if subset == "all":
            enzyme = default_enzyme
        else:
            raise ValueError("subset not configured: '{}'".format(subset))
    else:
        try:
            enzyme = subset_conf["enzyme"]
        except KeyError:
            enzyme = default_enzyme
    if enzyme not in config["software"]["percolator"]["enzymes"]:
        raise ValueError("unknown enzyme: '{}'".format(enzyme))
    return enzyme


def percolator_input_files(wc,config):
    group_to_samples = get_samples(wc.ds,wc.subset,wc.grouping,config)
    return ["out/{db}/split_pins/{ds}/{subset}/{sample}.{sdb}.pin".format(
            db=wc.db,ds=wc.ds,subset=wc.subset,sample=sample,sdb=wc.sdb)
            for sample in group_to_samples[wc.group]]


def sample_data_file(wildcards,config):
    ds = wildcards.ds
    sample = wildcards.sample
    ds_meta = dataset_metadata(ds,"all",config)
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
        ds_meta = get_dataset_metadata(ds,config)
        os.makedirs(ds_dir,exist_ok=True)
        with open(ds_meta_file,"w") as f:
            json.dump(ds_meta,f)


def sync_sample_proxy_files(ds,config):
    ds_dir = dataset_dir(ds,config)
    if os.path.exists(ds_dir):
        for item_name in os.listdir(ds_dir):
            item_path = os.path.join(ds_dir,item_name)
            if os.path.isfile(item_path) and item_path.endswith("_proxy.json"):
                os.remove(item_path)
    if dataset_source(ds,config) != "local":
        os.makedirs(ds_dir,exist_ok=True)
        ds_meta = dataset_metadata(ds,"all",config)
        for sample,sample_meta in ds_meta["samples"].items():
            file_url = sample_meta["url"]
            mod_dt = utc_strptime(sample_meta["last-modified"])
            file_mtime = timegm(mod_dt.timetuple())
            proxy_file_name = "{}_proxy.json".format(url_basename(file_url))
            proxy_file_path = os.path.join(ds_dir,proxy_file_name)
            with open(proxy_file_path,"w") as f:
                json.dump(sample_meta,f)
            os.utime(proxy_file_path,(time.time(),file_mtime))
