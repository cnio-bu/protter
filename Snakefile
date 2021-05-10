from functools import partial
import os
import re

import pandas as pd

from scripts.common import (get_samples,
                            load_sample_sheet,
                            pull_download_sheet)

from scripts.workflow import (comet_input_file,
                              dataset_groupings,
                              dataset_dbs,
                              dataset_subsets,
                              download_sample_file_checksum,
                              download_sample_file_url,
                              download_sample_output_pattern,
                              group_pin_files,
                              raw_input_file,
                              mzml_output_pattern,
                              percolator_enzyme)


configfile: os.path.join(workflow.basedir, "config.yaml")

config["dataset_path"] = os.path.relpath(os.path.join(
    workflow.basedir, config["dataset_path"]))
config["download_path"] = os.path.relpath(os.path.join(
    workflow.basedir, config["download_path"]))

sample_file = config["samples"]
samples = load_sample_sheet(sample_file)
downloads = pull_download_sheet(samples)

wildcard_constraints:
    dl_file = "|".join(re.escape(x) for x in downloads["dl_file"]),
    sample = "|".join(re.escape(x) for x in samples["sample"])

download_sample_file_checksum = partial(download_sample_file_checksum, downloads=downloads)
download_sample_file_url = partial(download_sample_file_url, downloads=downloads)
download_sample_output_pattern = partial(download_sample_output_pattern, config=config)
raw_input_file = partial(raw_input_file, config=config, samples=samples)
mzml_output_pattern = partial(mzml_output_pattern, config=config)
comet_input_file = partial(comet_input_file, config=config, samples=samples)
group_pin_files = partial(group_pin_files, samples=samples)
percolator_enzyme = partial(percolator_enzyme, config=config, samples=samples)


def input_crux_percolator():
    for ds in config["datasets"]:
        if not config["datasets"][ds]["enabled"]:
            continue
        if ds not in samples["dataset"].unique():
            continue
        for subset in dataset_subsets(ds,samples):
            for grouping in dataset_groupings(ds,config):
                for group in get_samples(ds,subset,grouping,samples):
                    for db in dataset_dbs(ds,config):
                        for em in config["software"]["percolator"]["modes"]:
                            for sdb in config["dbs"][db]["paths"]:
                                yield os.path.join("out",db,"percolator",ds,subset,em,grouping,
                                                   group,sdb,"percolator.target.psms.proc.txt")
    for ds in downloads["dataset"].unique():
        if not (ds in config["datasets"] and
                config["datasets"][ds]["enabled"]):
            continue
        yield os.path.join(config["download_path"],ds,"sha1checksums.txt")


rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        input_crux_percolator(),


include: "rules/download_sample.smk"
include: "rules/raw_to_mzml.smk"
include: "rules/procdb.smk"
include: "rules/comet.smk"
include: "rules/percolator.smk"
include: "rules/process_psms.smk"
