#!/usr/bin/env python

import os
from urllib.parse import unquote,urlparse


def dataset_source(ds,config):
    ds_conf = config["datasets"][ds]
    return ds_conf["source"] if "source" in ds_conf else "local"


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
