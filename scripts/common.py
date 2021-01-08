from calendar import timegm
from contextlib import contextmanager
from datetime import datetime,timezone
import functools
import gzip
import os
import re
from urllib.parse import unquote,urlparse,urlunparse

from dateutil.relativedelta import relativedelta
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
    ds_dir = dataset_dir(ds,config)

    ds_meta = {
        "config": config["datasets"][ds],
        "samples": {}
    }

    if ds_src == "PRIDE":

        pride_proj_meta = query_pride_project_metadata(ds,config)
        pride_file_meta = query_pride_file_metadata(ds,config)
        ds_url = get_pride_dataset_url(ds,pride_proj_meta,pride_file_meta)

        for rec in pride_file_meta["list"]:
            sample,file_ext,gzip_ext = split_gzip_ext(rec["fileName"])
            file_fmt = file_ext[1:]
            if ds_fmt_regex.match(file_fmt):
                if sample in ds_meta["samples"]:
                    raise ValueError("dataset '{}' has duplicate sample '{}'".format(ds,sample))
                ds_meta["samples"][sample] = {
                    "size": rec["fileSize"],
                    "url": ds_url + rec["fileName"]
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


def get_group_enzyme(ds,subset,grouping,group,samples,default=None):
    group_samples = get_group_sample_meta(ds,subset,grouping,group,samples)
    if "enzyme" in group_samples.columns:
        group_enzymes = set(group_samples["enzyme"].dropna())
        if len(group_enzymes) != 1:
            raise ValueError(
                "enzyme not configured correctly for group: '{}'".format(group))
        group_enzyme = group_enzymes.pop()
    else:
        group_enzyme = default
    return group_enzyme


def get_group_meta_value(ds,subset,grouping,group,samples,meta_key,default=None):
    group_samples = get_group_sample_meta(ds,subset,grouping,group,samples)
    if meta_key in group_samples.columns:
        meta_values = set(group_samples[meta_key].dropna())
        if len(meta_values) != 1:
            raise ValueError(
                "'{}' metadata not configured correctly for group: '{}'".format(meta_key,group))
        meta_value = meta_values.pop()
    else:
        meta_value = default
    return meta_value


def get_group_sample_meta(ds,subset,grouping,group,samples):
    group_samples = samples.xs(key=ds,level="dataset")
    group_to_samples = get_samples(ds,subset,grouping,samples)
    sample_ids = group_to_samples[group]
    group_samples = group_samples.loc[group_samples["sample"].isin(sample_ids)]
    return group_samples


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


def get_pride_dataset_url(ds,pride_proj_meta,pride_file_meta):

    # The dataset name is not guaranteed to be the same as the PRIDE
    # project accession, so we take that from the project metadata.
    proj_ac = pride_proj_meta["accession"]

    # Infer dataset URL from file metadata.
    ds_url = None
    for rec in pride_file_meta["list"]:
        file_url = rec["downloadLink"]
        parsed_url = urlparse(file_url)
        url_path_parts = parsed_url.path.split("/")
        pred_path = "/".join(url_path_parts[:-1]) + "/"
        pred_url_attrs = parsed_url[:2] + (pred_path,) + parsed_url[3:]
        pred_url = urlunparse(pred_url_attrs)
        if pred_url != ds_url:
            if ds_url is None:
                ds_url = pred_url
            else:
                raise ValueError(
                    "failed to infer PRIDE project URL of dataset: '{}'".format(ds))

    # We start from the assumption that the URL reflects the date of publication, but
    # allow for the possibility that it has been revised any time between then and now.
    pub_date = pride_proj_meta["publicationDate"]
    ds_date = datetime.strptime(pub_date,"%Y-%m-%d").date()
    date_now = datetime.utcnow().date()

    one_month = relativedelta(months=1)
    next_month = date_now + one_month
    next_month = next_month.replace(day=1)  # note to self: avoid future dates

    ds_url = None
    base_url = "http://ftp.pride.ebi.ac.uk/pride/data/archive"
    while ds_date < next_month:
        ds_url = "{}/{:d}/{:02d}/{}/".format(base_url,ds_date.year,ds_date.month,proj_ac)
        r = requests.head(ds_url)
        if r.status_code == 200:
            break
        ds_date = ds_date + one_month

    if ds_url is None:
        raise RuntimeError(
            "failed to determine PRIDE project URL of dataset '{}'".format(ds))

    return ds_url


def get_pride_file_mtime_info(file_urls):

    mtime_info = {}
    with requests.Session() as s:

        for file_url in file_urls:
            parsed_url = urlparse(file_url)

            if parsed_url.netloc != "ftp.pride.ebi.ac.uk":
                raise ValueError(
                    "invalid PRIDE file URL: '{}'".format(file_url))

            http_url = "http://{}{}".format(parsed_url.netloc,parsed_url.path)
            r = s.head(http_url)
            r.raise_for_status()

            mod_dt = datetime.strptime(r.headers["Last-Modified"],
                                       "%a, %d %b %Y %H:%M:%S GMT")
            mtime_info[file_url] = timegm(mod_dt.timetuple())

    return mtime_info


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


def load_config_file(config_file):
    with open(config_file) as f:
        return yaml.safe_load(f)


def load_sample_sheet(file):
    samples = pd.read_csv(file,sep="\t",dtype="string",na_values="NA",
                          keep_default_na=False)
    return samples.set_index(["dataset","sample"],drop=False,
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


def query_pride_file_metadata(ds,config):

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


def query_pride_project_metadata(ds,config):
    ds_conf = config["datasets"][ds]
    try:
        proj_ac = ds_conf["project_accession"]
    except KeyError:
        proj_ac = ds
    base_url = "http://www.ebi.ac.uk/pride/ws/archive"
    request_url = "{}/project/{}".format(base_url,proj_ac)
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
