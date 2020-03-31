#!/usr/bin/env python

from calendar import timegm
from datetime import datetime,timezone
from ftplib import FTP
import json
import os
import re
import time
from urllib.parse import urlparse,urlunparse

import requests

from .common import dataset_source,split_gzip_ext,url_basename


def _get_fmt_regex(fmts):
    pattern = "^({})$".format("|".join(fmts))
    return re.compile(pattern,re.IGNORECASE)


def _is_comet_fmt(fmt,config):
    regex = _get_fmt_regex(config["software"]["comet"]["fmts"])
    return regex.match(fmt) is not None


def _is_msconvert_fmt(fmt,config):
    regex = _get_fmt_regex(config["software"]["msconvert"]["fmts"])
    return regex.match(fmt) is not None


def _iter_list_file(file_path):
    with open(file_path) as f:
        for line in f:
            contents = line.strip()
            if contents:
                yield contents


def _make_ds_meta_file(ds,config,ds_meta_file):
    ds_fmt_regex = _get_fmt_regex([config["datasets"][ds]["fmt"]])
    ds_src = dataset_source(ds,config)
    ds_dir = dataset_dir(ds,config)

    ds_meta = {
        "config": config["datasets"][ds],
        "samples": {}
    }

    if ds_src == "PRIDE":

        file_meta = _pride_file_metadata(ds,config)

        for rec in file_meta["list"]:
            sample,file_ext,gzip_ext = split_gzip_ext(rec["fileName"])
            file_fmt = file_ext[1:]
            if ds_fmt_regex.match(file_fmt):
                assert sample not in ds_meta["samples"], \
                    "sample names must be unique within a PRIDE dataset"
                file_url = rec["downloadLink"]
                ds_meta["samples"][sample] = {
                    "size": rec["fileSize"],
                    "url": file_url
                }

        mtime_info = _pride_file_mtime_info(
            (x["url"] for x in ds_meta["samples"].values())
        )
        for sample,sample_meta in ds_meta["samples"].items():
            url = sample_meta["url"]
            mod_dt = datetime.fromtimestamp(mtime_info[url],timezone.utc)
            sample_meta["last-modified"] = _utc_strftime(mod_dt)

    elif ds_src == "local":

        if os.path.isdir(ds_dir):
            for item_name in sorted(os.listdir(ds_dir)):
                item_path = os.path.join(ds_dir,item_name)
                if os.path.isfile(item_path):
                    sample,file_ext,gzip_ext = split_gzip_ext(item_name)
                    file_fmt = file_ext[1:]
                    if (ds_fmt_regex.match(file_fmt) and
                            sample not in ds_meta["samples"]):
                        ds_meta["samples"][sample] = {
                            "path": item_path
                        }

    else:
        raise ValueError(
            "unsupported proteomics data source: '{}'".format(ds_src))

    os.makedirs(ds_dir,exist_ok=True)
    with open(ds_meta_file,"w") as f:
        json.dump(ds_meta,f)

    return ds_meta


def _pride_dataset_readme_url(ds,ds_meta):
    readme_url = None
    for sample_meta in ds_meta["samples"].values():
        file_url = sample_meta["url"]
        parsed_url = urlparse(file_url)
        url_path_parts = parsed_url.path.split("/")
        ds_dir_path = "/".join(url_path_parts[:-1]) + "/"
        pred_path = ds_dir_path + "README.txt"
        pred_url_attrs = parsed_url[:2] + (pred_path,) + parsed_url[3:]
        pred_url = urlunparse(pred_url_attrs)
        if pred_url != readme_url:
            if readme_url is None:
                readme_url = pred_url
            else:
                raise ValueError(
                    "cannot infer README URL for dataset: '{}'".format(ds))
    return readme_url


def _pride_file_metadata(ds,config):
    ds_conf = config["datasets"][ds]
    try:
        proj_ac = ds_conf["project_accession"]
    except KeyError:
        proj_ac = ds

    base_url = "http://www.ebi.ac.uk/pride/ws/archive"
    request_url = "{}/file/list/project/{}".format(base_url,proj_ac)
    response = requests.get(request_url)

    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        if response.status_code == 401:
            if "project_accession" in ds_conf:
                raise ValueError("dataset '{}' has invalid PRIDE project "
                                 "accession: '{}'".format(ds,proj_ac))
            else:
                raise ValueError("please specify a PRIDE project accession "
                                 "for dataset '{}'".format(ds))
        else:
            raise e

    return response.json()


def _pride_file_mtime_info(file_urls):
    pride_domain = "ftp.pride.ebi.ac.uk"
    mtime_info = {}
    with FTP(pride_domain) as ftp:
        ftp.login()

        for file_url in file_urls:

            parsed_url = urlparse(file_url)
            if parsed_url.netloc != pride_domain:
                raise ValueError(
                    "invalid PRIDE file URL: '{}'".format(file_url))

            response = ftp.voidcmd("MDTM {}".format(parsed_url.path))
            _,mdtm_ts = response.split()
            mod_dt = datetime.strptime(mdtm_ts, "%Y%m%d%H%M%S")
            mtime_info[file_url] = timegm(mod_dt.timetuple())

    return mtime_info


def _utc_strftime(dt):
    """Format datetime object as UTC timestamp string."""
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def _utc_strptime(ts):
    """Make datetime object from UTC timestamp string."""
    return datetime.strptime(ts,"%Y-%m-%dT%H:%M:%SZ")


def comet_input_file(wildcards,config):
    ds = wildcards.ds
    ds_fmt = config["datasets"][ds]["fmt"]

    if _is_comet_fmt(ds_fmt,config):
        input_file = sample_data_file(wildcards,config)
    elif _is_msconvert_fmt(ds_fmt,config):
        input_file = msconvert_output_pattern(config)
    else:
        raise ValueError(
            "unsupported proteomics file format: '{}'".format(ds_fmt))
    return input_file


def dataset_dbs(ds,config):
    return [db for db in config["dbs"]
            if config["dbs"][db]["enabled"]
            and db in config["datasets"][ds]["dbs"]]


def dataset_dir(ds,config):
    return os.path.join(config["dataset_path"],ds)


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

    if _is_msconvert_fmt(ds_fmt,config):
        msconvert_input = sample_data_file(wildcards,config)
    elif _is_comet_fmt(ds_fmt,config):  # i.e. no need to run msconvert
        msconvert_input = ""
    else:
        raise ValueError(
            "unsupported proteomics file format: '{}'".format(ds_fmt))

    return msconvert_input


def msconvert_output_pattern(config):
    return os.path.join(config["dataset_path"],"{ds}","{sample}.mzML.gz")


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
                readme_url = _pride_dataset_readme_url(ds,ds_meta)
                mtime_info = _pride_file_mtime_info([readme_url])
                latest_mtime = mtime_info[readme_url]
            else:
                raise ValueError(
                    "unsupported proteomics data source: '{}'".format(ds_src))

            update_needed = meta_file_mtime < latest_mtime

    if update_needed:
        _make_ds_meta_file(ds,config,ds_meta_file)


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
            mod_dt = _utc_strptime(sample_meta["last-modified"])
            file_mtime = timegm(mod_dt.timetuple())
            proxy_file_name = "{}_proxy.json".format(url_basename(file_url))
            proxy_file_path = os.path.join(ds_dir,proxy_file_name)
            with open(proxy_file_path,"w") as f:
                json.dump(sample_meta,f)
            os.utime(proxy_file_path,(time.time(),file_mtime))
