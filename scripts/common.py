import os
import re
from urllib.parse import unquote,urlparse

import pandas as pd
from ratelimiter import RateLimiter
from requests.adapters import HTTPAdapter
from requests.exceptions import HTTPError
from requests import Session
import yaml


def _get_fmt_regex(fmts):
    pattern = "^({})$".format("|".join(re.escape(x) for x in fmts))
    return re.compile(pattern,re.IGNORECASE)


def dataset_dir(ds,config):
    return os.path.join(config["dataset_path"],ds)


def dataset_source(ds,config):
    ds_conf = config["datasets"][ds]
    return ds_conf["source"] if "source" in ds_conf else "local"


def download_dir(ds,config):
    return os.path.join(config["download_path"],ds)


def get_dataset_metadata(ds,config):

    ds_fmt_regex = _get_fmt_regex([config["datasets"][ds]["fmt"]])
    ds_src = dataset_source(ds,config)

    ds_meta = dict()
    if ds_src == "PRIDE":

        if "projects" in config["datasets"][ds]:
            proj_accs = config["datasets"][ds]["projects"]
        elif "project" in config["datasets"][ds]:
            proj_accs = [config["datasets"][ds]["project"]]
        else:
            raise ValueError(
                "no PRIDE project specified for dataset: '{}'".format(ds))

        with Session() as session:
            session.mount("https://www.ebi.ac.uk/pride/",
                          HTTPAdapter(max_retries=3))

            for proj_acc in proj_accs:
                proj_file_meta = query_pride_file_metadata(session,proj_acc)

                for rec in proj_file_meta:

                    if proj_acc not in rec["projectAccessions"]:
                        raise ValueError(
                            "PRIDE file metadata lacks project accession: '{}'".format(proj_acc))

                    sample,file_ext,gzip_ext = split_gzip_ext(rec["fileName"])
                    file_fmt = file_ext[1:]
                    if ds_fmt_regex.match(file_fmt):
                        if sample in ds_meta:
                            raise ValueError(
                                "dataset '{}' has duplicate sample '{}'".format(ds,sample))

                        file_url = None
                        for file_loc_meta in rec["publicFileLocations"]:
                            if file_loc_meta["name"] == "FTP Protocol":
                                file_url = file_loc_meta["value"]
                                break

                        if file_url is None:
                            raise ValueError(
                                "FTP URL not found for sample '{}'".format(sample))

                        ds_meta[sample] = {
                            "project": proj_acc,
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
                            sample not in ds_meta):
                        ds_meta[sample] = {
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

    if "subset" in rel_samples.columns:
        if subset == "all" and set(rel_samples["subset"]) != {"all"}:
            raise ValueError(
                "subset 'all' not configured correctly for dataset '{}'".format(ds))
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


def is_raw_fmt(fmt,config):
    regex = _get_fmt_regex(["raw"])
    return regex.match(fmt) is not None


def is_remote_url(file):
    remote_url_prefix = re.compile("(?:ftp|http|https)://",
                                   re.IGNORECASE)
    return remote_url_prefix.match(file) is not None


def load_config_file(config_file):
    with open(config_file) as f:
        return yaml.safe_load(f)


def load_sample_sheet(file):
    samples = pd.read_csv(file,sep="\t",dtype="string",na_values="NA",
                          keep_default_na=False)
    return samples.set_index(["dataset","sample"],drop=False,
                             verify_integrity=True)


def pull_download_sheet(samples):

    downloads = samples.loc[samples["file"].apply(is_remote_url),:]

    if len(downloads) > 0:

        if "subset" in samples.columns:
            downloads = downloads.loc[~downloads["subset"].isna(),:]

        if downloads["checksum"].isna().any():
            raise ValueError("all download files must have an associated checksum")

        downloads = downloads[["dataset","file","checksum"]].assign(
            dl_file=downloads["file"].apply(url_basename))

        downloads = downloads.set_index(["dataset","dl_file"],drop=False,
                                        verify_integrity=True)
    else:
        downloads = pd.DataFrame(columns=["dataset","file","checksum","dl_file"
                                          ]).set_index(["dataset","dl_file"], drop=False)

    return downloads


@RateLimiter(max_calls=3, period=1)
def query_pride_file_metadata(session,proj_acc):

    base_url = "https://www.ebi.ac.uk/pride/ws/archive/v2"
    request_url = "{}/files/byProject?accession={}".format(base_url,proj_acc)
    response = session.get(request_url)

    try:
        response.raise_for_status()
    except HTTPError as e:
        if response.status_code == 401:
            raise ValueError("invalid PRIDE project accession: '{}'".format(proj_acc))
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
