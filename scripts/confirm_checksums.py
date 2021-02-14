from contextlib import redirect_stdout,redirect_stderr
import os

from common import get_dataset_metadata,load_config_file,split_gzip_ext


input_cksum_files = snakemake.input.cksum_files
config_file = snakemake.input.config_file
output_cksum_file = snakemake.output.cksum_file
out_log_file = snakemake.log.o
err_log_file = snakemake.log.e
ds = snakemake.wildcards.ds


with open(out_log_file,"w") as out_log_f, open(err_log_file,"w") as err_log_f:
    with redirect_stdout(out_log_f), redirect_stderr(err_log_f):
        config = load_config_file(config_file)
        ds_meta = get_dataset_metadata(ds,config)

        with open(output_cksum_file,"w") as out_f:
            for input_cksum_file in input_cksum_files:

                dl_file_path,input_cksum_ext = os.path.splitext(input_cksum_file)
                if input_cksum_ext != ".sha1":
                    raise ValueError(
                        "checksum file has unknown extension: '{}'".format(input_cksum_ext))

                dl_file_name = os.path.basename(dl_file_path)
                sample,file_ext,gzip_ext = split_gzip_ext(dl_file_name)

                exp_cksum = ds_meta["samples"][sample]["checksum"]
                with open(input_cksum_file) as in_f:
                    obs_cksum = in_f.read()
                if obs_cksum != exp_cksum:
                    raise ValueError(
                        "checksum mismatch for file: '{}'".format(dl_file_path))

                out_f.write("{}  {}\n".format(obs_cksum,dl_file_name))
