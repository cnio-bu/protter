from contextlib import redirect_stdout,redirect_stderr
import os

from common import get_dataset_metadata,load_config_file, \
                   split_gzip_ext,url_basename

input_cksum_files = snakemake.input.cksum_files
config_file = snakemake.input.config_file
output_cksum_file = snakemake.output.cksum_file
out_log_file = snakemake.log.o
err_log_file = snakemake.log.e
ds = snakemake.wildcards.ds

dl_dirs = set()
input_cksum_file_names = list()
for cksum_file_path in input_cksum_files:
    cksum_file_dir,cksum_file_name = os.path.split(cksum_file_path)
    input_cksum_file_names.append(cksum_file_name)
    dl_dirs.add(cksum_file_dir)

if len(dl_dirs) != 1:
    raise ValueError("inconsistent download directory")
dl_dir = dl_dirs.pop()

# We want to confirm checksums for all files
# that have been downloaded for this dataset.
obs_cksum_file_names = set([x for x in os.listdir(dl_dir)
                            if x.endswith(".sha1")])

with open(out_log_file,"w") as out_log_f, open(err_log_file,"w") as err_log_f:
    with redirect_stdout(out_log_f), redirect_stderr(err_log_f):
        config = load_config_file(config_file)
        ds_meta = get_dataset_metadata(ds,config)

        exp_cksum_file_names = set()
        with open(output_cksum_file,"w") as out_f:
            for sample in sorted(ds_meta.keys()):
                sample_meta = ds_meta[sample]

                sample_file_name = url_basename(sample_meta["file"])
                cksum_file_name = "{}.sha1".format(sample_file_name)
                exp_cksum_file_names.add(cksum_file_name)

                if cksum_file_name in obs_cksum_file_names:
                    cksum_file_path = os.path.join(dl_dir,cksum_file_name)
                    dl_file_path = os.path.join(dl_dir,sample_file_name)

                    exp_cksum = ds_meta[sample]["checksum"]
                    with open(cksum_file_path) as in_f:
                        obs_cksum = in_f.read()
                    if obs_cksum != exp_cksum:
                        raise ValueError(
                            "checksum mismatch for file: '{}'".format(dl_file_path))

                    out_f.write("{}  {}\n".format(obs_cksum,sample_file_name))

        for cksum_file_name in input_cksum_file_names:
            if cksum_file_name not in exp_cksum_file_names:
                raise ValueError(
                    "unexpected checksum file: '{}'".format(cksum_file_name))
