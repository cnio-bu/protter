#!/usr/bin/env python3
"""Gather Percolator PSM files."""

from argparse import ArgumentParser
from collections import defaultdict
from glob import iglob
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from common import load_config_file


if __name__ == "__main__":

    ap = ArgumentParser(description=__doc__)
    ap.add_argument("config_file",
                    help="Input config YAML file indicating the database(s) and"
                         " dataset(s) for which PSM files should be gathered.")
    ap.add_argument("output_file",
                    help="Output ZIP file containing the gathered PSM files.")
    args = ap.parse_args()

    config_file = args.config_file
    output_file = args.output_file


    config = load_config_file(config_file)

    db_to_datasets = defaultdict(list)
    for ds in config["datasets"]:
        if not config["datasets"][ds].get("enabled",False):
            continue
        for db in config["datasets"][ds]["dbs"]:
            if not config["dbs"][db].get("enabled",False):
                continue
            db_to_datasets[db].append(ds)

    script_dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.dirname(script_dir)
    out_dir = os.path.join(base_dir,"out")

    with TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)

        for db,datasets in db_to_datasets.items():
            for ds in datasets:
                print("Copying PSM files for searches of dataset"
                      " '{}' against database '{}'...".format(ds,db))
                ds_dir = tmp_dir / db / ds
                os.makedirs(ds_dir,exist_ok=True)

                search_patt = "{}/{}/percolator/{}/**/percolator.target.psms.proc.txt".format(
                    out_dir,db,ds)
                for in_file_glob_path in iglob(search_patt,recursive=True):
                    in_file_rel_path = Path(os.path.relpath(in_file_glob_path,start=out_dir))
                    db,_,ds,subset,em,grouping,group,sdb,file_name = in_file_rel_path.parts
                    out_file_name = "{}.{}.{}.{}.{}.percolator.psms.proc.txt".format(subset,em,
                                                                                     grouping,
                                                                                     group,sdb)
                    out_file_path = ds_dir / out_file_name
                    shutil.copyfile(in_file_glob_path,out_file_path)

        print("Making output archive file...")
        archive_base,archive_ext = os.path.splitext(output_file)
        archive_fmt = archive_ext[1:]
        shutil.make_archive(archive_base,archive_fmt,tmp_dir)

    print("Done.")
