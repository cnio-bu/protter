#!/usr/bin/env python

import os
from urllib.parse import urlparse,urlunparse

import requests
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from .common import dataset_source,split_gzip_ext,url_basename


def _get_input_file_by_format(fmt,wildcards,config):
    ds = wildcards.ds
    sample = wildcards.sample
    ds_src = dataset_source(ds,config)
    exp_file_name = "{}.{}".format(sample,fmt.lower())
    exp_gzip_name = "{}.gz".format(exp_file_name)
    basename = os.path.basename if ds_src == "local" else url_basename
    for file in dataset_input_files(ds,config,sort=True):
        file_stem,file_ext,gzip_ext = split_gzip_ext(basename(file))
        file_name = file_stem + file_ext.lower() + gzip_ext.lower()
        if file_name in (exp_file_name,exp_gzip_name):
            input_file = file
            break
    if ds_src != "local":
        input_file = _get_remote_input(input_file,ds_src)
    return input_file


def _get_lc_file_format(file_path):
    file_ext = split_gzip_ext(file_path)[1]
    if not file_ext:
        raise ValueError("cannot infer format of file: '{}'".format(file_path))
    return file_ext[1:].lower()


def _get_remote_input(file_pattern,source):
    if source == "PRIDE":
        # Though PRIDE files are hosted on an FTP site, they
        # are easier to access via the HTTP remote provider.
        provider = HTTPRemoteProvider()
        remote_patt = provider.remote(file_pattern,insecure=True)
    else:
        raise ValueError(
            "unsupported proteomics data source: '{}'".format(source))
    return remote_patt


def _strip_scheme_prefix(uri):
    # Stripping leading slashes is a hack, but a necessary one,
    # as urlunparse appears to return a protocol-relative URL.
    return urlunparse(("",) + urlparse(uri)[1:]).lstrip("/")


def comet_input_file(wildcards,config):
    ds = wildcards.ds
    sample = wildcards.sample
    ds_fmt = config["datasets"][ds]["fmt"]
    lc_ds_fmt = ds_fmt.lower()

    comet_fmts = [x.lower()
        for x in config["software"]["comet"]["fmts"]]
    msconvert_fmts = [x.lower()
        for x in config["software"]["msconvert"]["fmts"]]

    if lc_ds_fmt in comet_fmts:
        input_file = _get_input_file_by_format(ds_fmt,wildcards,config)
    elif lc_ds_fmt in msconvert_fmts:
        input_file = msconvert_output_pattern(config)
    else:
        raise ValueError(
            "unsupported proteomics file format: '{}'".format(ds_fmt))
    return input_file


def dataset_input_files(ds,config,sort=False):
    """Get workflow input files for the specified dataset.

    Returns file paths if the dataset source is local,
    otherwise returns file URLs.
    """
    ds_conf = config["datasets"][ds]
    ds_src = dataset_source(ds,config)
    lc_ds_fmt = ds_conf["fmt"].lower()

    input_files = []
    if ds_src == "PRIDE":

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

        response_data = response.json()
        for rec in response_data["list"]:
            file_name = rec["fileName"]
            lc_file_fmt = _get_lc_file_format(rec["fileName"])
            if lc_file_fmt == lc_ds_fmt:
                file_uri = _strip_scheme_prefix(rec["downloadLink"])
                input_files.append(file_uri)

    elif ds_src == "local":

        ds_path = os.path.join(config["dataset_path"],ds)
        if os.path.isdir(ds_path):
            for file_name in os.listdir(ds_path):
                file_path = os.path.join(ds_path,file_name)
                if os.path.isfile(file_path):
                    lc_file_fmt = _get_lc_file_format(file_path)
                    if lc_file_fmt == lc_ds_fmt:
                        input_files.append(file_path)

    else:
        raise ValueError(
            "unsupported proteomics data source: '{}'".format(ds_src))

    return sorted(input_files) if sort else input_files


def msconvert_input_file(wildcards,config):
    ds = wildcards.ds
    sample = wildcards.sample
    ds_fmt = config["datasets"][ds]["fmt"]
    lc_ds_fmt = ds_fmt.lower()

    comet_fmts = [x.lower()
        for x in config["software"]["comet"]["fmts"]]
    msconvert_fmts = [x.lower()
        for x in config["software"]["msconvert"]["fmts"]]

    if lc_ds_fmt in msconvert_fmts:
        msconvert_input = _get_input_file_by_format(ds_fmt,wildcards,config)
    elif lc_ds_fmt in comet_fmts:  # i.e. no need to run msconvert
        msconvert_input = ""
    else:
        raise ValueError(
            "unsupported proteomics file format: '{}'".format(ds_fmt))

    return msconvert_input


def msconvert_output_pattern(config):
    return os.path.join(config["dataset_path"],"{ds}","{sample}.mzML.gz")
