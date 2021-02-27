from contextlib import redirect_stdout,redirect_stderr
import gzip


gzip_pin_files = snakemake.input.pins
text_pin_file = snakemake.output.pin
out_log_file = snakemake.log.o
err_log_file = snakemake.log.e

with open(out_log_file,"w") as out_log_f, open(err_log_file,"w") as err_log_f:
    with redirect_stdout(out_log_f), redirect_stderr(err_log_f):
        exp_hdr_line = None
        with open(text_pin_file,"wt") as out_f:
            for gzip_pin_file in gzip_pin_files:
                with gzip.open(gzip_pin_file,"rt") as in_f:

                    hdr_line = next(in_f)
                    if hdr_line != exp_hdr_line:
                        if exp_hdr_line is None:
                            exp_hdr_line = hdr_line
                            out_f.write(hdr_line)
                        else:
                            raise ValueError("PIN header line mismatch")

                    for line in in_f:
                        out_f.write(line)
