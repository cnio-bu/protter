from contextlib import contextmanager
import functools
import gzip
import os
import re
from urllib.parse import unquote,urlparse

import pandas as pd
import requests
import yaml


def _get_fmt_regex(fmts):
    pattern = "^({})$".format("|".join(fmts))
    return re.compile(pattern,re.IGNORECASE)


def dataset_dir(ds,config):
    return os.path.join(config["dataset_path"],ds)


def dataset_source(ds,config):
    ds_conf = config["datasets"][ds]
    return ds_conf["source"] if "source" in ds_conf else "local"


def get_dataset_metadata(ds,config):

    ds_fmt_regex = _get_fmt_regex([config["datasets"][ds]["fmt"]])
    ds_src = dataset_source(ds,config)

    ds_meta = {
        "samples": {}
    }

    if ds_src == "PRIDE":

        proj_acn = config["datasets"][ds]["project_accession"]
        pride_file_meta = query_pride_file_metadata(proj_acn)

        for rec in pride_file_meta:
            sample,file_ext,gzip_ext = split_gzip_ext(rec["fileName"])
            file_fmt = file_ext[1:]
            if ds_fmt_regex.match(file_fmt):
                if sample in ds_meta["samples"]:
                    raise ValueError("PRIDE project '{}' has duplicate sample '{}'".format(proj_acn,sample))

                file_url = None
                for file_loc_meta in rec["publicFileLocations"]:
                    if file_loc_meta["name"] == "FTP Protocol":
                        file_url = file_loc_meta["value"]
                        break

                if file_url is None:
                    raise ValueError("FTP URL not found for sample '{}'".format(sample))

                ds_meta["samples"][sample] = {
                    "size": rec["fileSizeBytes"],
                    "checksum": rec["checksum"],
                    "file": file_url
                }

    elif ds_src == "local":

        ds_dir = dataset_dir(ds,config)
        if os.path.isdir(ds_dir):
            for item_name in sorted(os.listdir(ds_dir)):
                item_path = os.path.join(ds_dir,item_name)
                if os.path.isfile(item_path):
                    sample,file_ext,gzip_ext = split_gzip_ext(item_name)
                    file_fmt = file_ext[1:]
                    if (ds_fmt_regex.match(file_fmt) and
                            sample not in ds_meta["samples"]):
                        ds_meta["samples"][sample] = {
                            "file": item_path
                        }
    else:
        raise ValueError(
            "unsupported proteomics data source: '{}'".format(ds_src))

    return ds_meta


def get_group_enzyme(ds,subset,grouping,group,samples,default=None):
    group_samples = get_group_sample_meta(ds,subset,grouping,group,samples)

    group_enzyme = default
    if "enzyme" in group_samples.columns:
        group_enzymes = set(group_samples["enzyme"].dropna())
        if len(group_enzymes) == 1:
            group_enzyme = group_enzymes.pop()
        elif len(group_enzymes) > 1:
            raise ValueError(
                "enzyme not configured correctly for group: '{}'".format(group))

    return group_enzyme


def get_group_meta_value(ds,subset,grouping,group,samples,meta_key,default=None):
    group_samples = get_group_sample_meta(ds,subset,grouping,group,samples)

    meta_value = default
    if meta_key in group_samples.columns:
        meta_values = set(group_samples[meta_key].dropna())
        if len(meta_values) == 1:
            meta_value = meta_values.pop()
        elif len(meta_values) > 1:
            raise ValueError(
                "'{}' metadata not configured correctly for group: '{}'".format(meta_key,group))

    return meta_value


def get_group_sample_meta(ds,subset,grouping,group,samples):
    group_samples = samples.xs(key=ds,level="dataset")
    group_to_samples = get_samples(ds,subset,grouping,samples)
    sample_ids = group_to_samples[group]
    group_samples = group_samples.loc[group_samples["sample"].isin(sample_ids)]
    return group_samples


def get_samples(ds,subset,grouping,samples):
    rel_samples = samples.xs(key=ds,level="dataset")

    if subset != "all":
        if "subset" not in rel_samples.columns:
            raise ValueError("subset not configured: '{}'".format(subset))
        rel_samples = rel_samples.loc[rel_samples["subset"] == subset]

    group_to_samples = dict()
    if grouping in rel_samples.columns:
        for group,sample in zip(rel_samples[grouping],
                                rel_samples["sample"]):
            try:
                group_to_samples[group].append(sample)
            except KeyError:
                group_to_samples[group] = [sample]
    elif grouping == "single":  # i.e. group == sample
        for sample in rel_samples["sample"]:
            group_to_samples[sample] = [sample]
    else:
        raise ValueError("grouping not configured: '{}'".format(grouping))

    return group_to_samples


def is_comet_fmt(fmt,config):
    regex = _get_fmt_regex(config["software"]["comet"]["fmts"])
    return regex.match(fmt) is not None


def is_msconvert_fmt(fmt,config):
    regex = _get_fmt_regex(config["software"]["msconvert"]["fmts"])
    return regex.match(fmt) is not None


def is_wget_url(file):
    wget_url_prefix = re.compile("(?:ftp|http|https)://",
                                 re.IGNORECASE)
    return wget_url_prefix.match(file) is not None


def load_config_file(config_file):
    with open(config_file) as f:
        return yaml.safe_load(f)


def load_sample_sheet(file):
    samples = pd.read_csv(file,sep="\t",dtype="string",na_values="NA",
                          keep_default_na=False)
    return samples.set_index(["dataset","sample"],drop=False,
                             verify_integrity=True)


def pull_download_sheet(samples):

    downloads = samples.loc[samples["file"].apply(is_wget_url),:]

    if downloads["checksum"].isna().any():
        raise ValueError("all download files must have an associated checksum")

    downloads = downloads[["dataset","file","checksum"]].assign(
        dl_file=downloads["file"].apply(url_basename))

    return downloads.set_index(["dataset","dl_file"],drop=False,
                               verify_integrity=True)


@contextmanager
def open_as_text(file,mode="r"):
    if mode not in {"a","r","w","x"}:
        raise ValueError("unknown file mode: '{}'".format(mode))
    text_mode = mode + "t"
    if str(file).lower().endswith(".gz"):
        opener = functools.partial(gzip.open,mode=text_mode)
    else:
        opener = functools.partial(open,mode=text_mode)
    file_obj = None
    try:
        file_obj = opener(file)
        yield file_obj
    finally:
        if file_obj is not None:
            file_obj.close()


def query_pride_file_metadata(proj_acn):

    base_url = "https://www.ebi.ac.uk/pride/ws/archive/v2"
    request_url = "{}/files/byProject?accession={}".format(base_url,proj_acn)
    response = requests.get(request_url)

    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        if response.status_code == 401:
            raise ValueError("invalid PRIDE project accession: '{}'".format(proj_acn))
        else:
            raise e

    return response.json()


def split_gzip_ext(file_path):
    parts = os.path.splitext(file_path)
    if parts[-1].lower() == ".gz":
        parts = os.path.splitext(parts[0]) + parts[-1:]
    else:
        parts += ("",)
    return parts


def url_basename(url):
    parsed_url = urlparse(url,allow_fragments=False)
    parsed_path = unquote(parsed_url.path)
    url_path_parts = parsed_path.split("/")
    basename = url_path_parts[-1]
    if not basename:
        raise ValueError("failed to get basename of URL: '{}'".format(url))
    return basename
