from functools import partial
import os

from scripts.common import dataset_source, split_gzip_ext, url_basename
from scripts.workflow import (comet_input_file, dataset_input_files,
                              msconvert_input_file, msconvert_output_pattern)


configfile: os.path.join(workflow.basedir, "config.yaml")

config["dataset_path"] = os.path.relpath( os.path.join(
    workflow.basedir, config["dataset_path"]))

msconvert_input_file = partial(msconvert_input_file, config=config)
msconvert_output_pattern = partial(msconvert_output_pattern, config=config)
comet_input_file = partial(comet_input_file, config=config)


def get_samples(ds,grouping):
    fpath = "res/data/prot/{ds}".format(ds=ds)
    #get the grouping statement from the config file and create a "gf" function with it
    exec('gf = lambda x: {}'.format(config["datasets"][ds]["groupings"][grouping]), globals())
    groups = {}
    ds_src = dataset_source(ds, config)
    basename = os.path.basename if ds_src == "local" else url_basename
    for fname in dataset_input_files(ds,config):
        f = split_gzip_ext(basename(fname))[0]
        g = gf(f)
        try:
            groups[g].append(f)
        except KeyError:
            groups[g] = [f]
    return groups

def input_crux_percolator():
    for ds in config["datasets"]:
        if config["datasets"][ds]["enabled"]:
            if "groupings" not in config["datasets"][ds]:
                config["datasets"][ds]["groupings"] = config["grouping_default"]
            for grouping in config["datasets"][ds]["groupings"]:
                for group in get_samples(ds,grouping):
                    for db in config["dbs"]:
                        if config["dbs"][db]["enabled"] and db in config["datasets"][ds]["dbs"]:
                            for em in config["software"]["percolator"]["modes"]:
                                for sdb in config["dbs"][db]["paths"]:
                                    yield "out/{db}/percolator/{ds}/{em}/{grouping}/{group}/{sdb}/percolator.target.peptides.txt".format(grouping=grouping,group=group,db=db,ds=ds,em=em,sdb=sdb)

rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        input_crux_percolator(),

include: "rules/msconvert.smk"
include: "rules/procdb.smk"
include: "rules/comet.smk"
include: "rules/percolator.smk"
