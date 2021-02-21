from contextlib import redirect_stdout,redirect_stderr
import gzip
import shutil


gzip_pin_file = snakemake.input.pin
text_pin_file = snakemake.output.pin
out_log_file = snakemake.log.o
err_log_file = snakemake.log.e

with open(out_log_file,"w") as out_log_f, open(err_log_file,"w") as err_log_f:
    with redirect_stdout(out_log_f), redirect_stderr(err_log_f):
        with gzip.open(gzip_pin_file,"rt") as in_f, open(text_pin_file,"wt") as out_f:
            shutil.copyfileobj(in_f,out_f)
