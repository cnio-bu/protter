import os,glob

configfile: "config.yaml"

def get_samples(ds,grouping):
    fpath = "res/data/prot/{ds}".format(ds=ds)
    #get the grouping statement from the config file and create a "gf" function with it
    exec('gf = lambda x: {}'.format(config["datasets"][ds]["groupings"][grouping]), globals())
    groups = {}
    for fname in glob.glob("{}/*.mzML".format(fpath)):
        f = os.path.splitext(os.path.basename(fname))[0]
        g = gf(f)
        try:
            groups[g].append(f)
        except KeyError:
            groups[g] = [f]
    return groups

def input_crux_percolator():
    for software in config["software"]["search"]:
        if config["software"]["search"][software]["enabled"]:
            for ds in config["datasets"]:
                if config["datasets"][ds]["enabled"]:
                    if "groupings" not in config["datasets"][ds]:
                        config["datasets"][ds]["groupings"] = config["grouping_default"]
                    for grouping in config["datasets"][ds]["groupings"]:
                        for group in get_samples(ds,grouping):
                            for db in config["dbs"]:
                                if config["dbs"][db]["enabled"] and db in config["datasets"][ds]["dbs"]:
                                    for em in ["tdc","mix-max"]:
                                        yield "out/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}/percolator.target.peptides.txt".format(grouping=grouping,group=group,sw=software,db=db,ds=ds,em=em)
'''
def input_add_confidence():
    for software in config["software"]["search"]:
        if config["software"]["search"][software]["enabled"]:
            for ds in config["datasets"]:
                if config["datasets"][ds]["enabled"]:
                    for db in config["dbs"]:
                        if "add_confidence" in config["dbs"][db]:
                            if config["dbs"][db]["enabled"] and db in config["datasets"][ds]["dbs"]:
                                    for sample in get_samples(ds,'single'):
                                        yield "out/{db}/add_confidence/{sw}/{ds}/{sample}.tsv".format(sw=software,db=db,ds=ds,sample=sample)
'''
rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        input_crux_percolator(),
        #input_add_confidence()

include: "rules/procdb.smk"
include: "rules/comet.smk"
#include: "rules/addconfidence.smk"
#include: "rules/xtandem.smk"
include: "rules/percolator.smk"
