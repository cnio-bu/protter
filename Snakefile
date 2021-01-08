from functools import partial
import os

from scripts.common import get_samples,load_sample_sheet
from scripts.workflow import (comet_input_file,
                              dataset_groupings,
                              dataset_dbs,
                              dataset_subsets,
                              download_sample_output_pattern,
                              msconvert_input_file,
                              msconvert_output_pattern,
                              msconvert_rule_config,
                              percolator_enzyme,
                              percolator_input_files,
                              sync_dataset_metadata,
                              sync_sample_proxy_files)


configfile: os.path.join(workflow.basedir, "config.yaml")

config["dataset_path"] = os.path.relpath(os.path.join(
    workflow.basedir, config["dataset_path"]))

sample_file = config["samples"]
samples = load_sample_sheet(sample_file)

download_sample_output_pattern = partial(download_sample_output_pattern, config=config)
msconvert_input_file = partial(msconvert_input_file, config=config, samples=samples)
msconvert_output_pattern = partial(msconvert_output_pattern, config=config)
msconvert_rule_config = partial(msconvert_rule_config, config=config)
comet_input_file = partial(comet_input_file, config=config, samples=samples)
percolator_input_files = partial(percolator_input_files, samples=samples)
percolator_enzyme = partial(percolator_enzyme, config=config, samples=samples)


def input_crux_percolator():
    for ds in config["datasets"]:
        if not config["datasets"][ds]["enabled"]:
            continue
        sync_dataset_metadata(ds,config)
        sync_sample_proxy_files(ds,config,samples)
        for subset in dataset_subsets(ds,samples):
            for grouping in dataset_groupings(ds,config):
                for group in get_samples(ds,subset,grouping,samples):
                    for db in dataset_dbs(ds,config):
                        for em in config["software"]["percolator"]["modes"]:
                            for sdb in config["dbs"][db]["paths"]:
                                yield os.path.join("out",db,"percolator",ds,subset,em,grouping,
                                                   group,sdb,"percolator.target.psms.proc.txt")


rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        input_crux_percolator(),


include: "rules/download_sample.smk"
include: "rules/msconvert.smk"
include: "rules/procdb.smk"
include: "rules/comet.smk"
include: "rules/percolator.smk"
include: "rules/process_psms.smk"
