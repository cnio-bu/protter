from calendar import timegm
from datetime import datetime,timezone
import os
import re
from urllib.parse import unquote,urlparse,urlunparse

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
    ds_dir = dataset_dir(ds,config)

    ds_meta = {
        "config": config["datasets"][ds],
        "samples": {}
    }

    if ds_src == "PRIDE":

        file_meta = get_pride_file_metadata(ds,config)

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

        mtime_info = get_pride_file_mtime_info(
            (x["url"] for x in ds_meta["samples"].values())
        )
        for sample,sample_meta in ds_meta["samples"].items():
            url = sample_meta["url"]
            mod_dt = datetime.fromtimestamp(mtime_info[url],timezone.utc)
            sample_meta["last-modified"] = utc_strftime(mod_dt)

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

    return ds_meta


def get_pride_dataset_readme_url(ds,ds_meta):
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


def get_pride_file_metadata(ds,config):

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


def get_pride_file_mtime_info(file_urls):

    mtime_info = {}
    for file_url in file_urls:
        parsed_url = urlparse(file_url)

        if parsed_url.netloc != "ftp.pride.ebi.ac.uk":
            raise ValueError(
                "invalid PRIDE file URL: '{}'".format(file_url))

        http_url = "http://{}{}".format(parsed_url.netloc,parsed_url.path)
        r = requests.head(http_url)
        r.raise_for_status()

        mod_dt = datetime.strptime(r.headers["Last-Modified"],
                                   "%a, %d %b %Y %H:%M:%S GMT")
        mtime_info[file_url] = timegm(mod_dt.timetuple())

    return mtime_info


def is_comet_fmt(fmt,config):
    regex = _get_fmt_regex(config["software"]["comet"]["fmts"])
    return regex.match(fmt) is not None


def is_msconvert_fmt(fmt,config):
    regex = _get_fmt_regex(config["software"]["msconvert"]["fmts"])
    return regex.match(fmt) is not None


def load_config_file(config_file):
    with open(config_file) as f:
        return yaml.safe_load(f)


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


def utc_strftime(dt):
    """Format datetime object as UTC timestamp string."""
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def utc_strptime(ts):
    """Make datetime object from UTC timestamp string."""
    return datetime.strptime(ts,"%Y-%m-%dT%H:%M:%SZ")
