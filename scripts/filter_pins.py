from contextlib import redirect_stdout,redirect_stderr
import re

in_pin_file = snakemake.input.pin
out_pin_file = snakemake.output.pin
out_log_file = snakemake.log.o
err_log_file = snakemake.log.e
peptide_length_range = snakemake.params.peptide_length_range


with open(out_log_file,"w") as out_log_f, open(err_log_file,"w") as err_log_f:
    with redirect_stdout(out_log_f), redirect_stderr(err_log_f):

        pep_len_rng_regex = re.compile("^(?P<min_pep_len>[0-9]+)\s+(?P<max_pep_len>[0-9]+)$")
        match = pep_len_rng_regex.match(peptide_length_range)
        groups = match.groupdict()
        min_pep_len = int(groups["min_pep_len"])
        max_pep_len = int(groups["max_pep_len"])

        with open(in_pin_file,"r") as in_f, open(out_pin_file,"w") as out_f:

            header_line = next(in_f)
            out_f.write(header_line)

            field_names = header_line.rstrip().split("\t")
            pep_len_idx = field_names.index("PepLen")

            for line in in_f:
                fields = line.rstrip().split("\t")
                pep_len = int(fields[pep_len_idx])
                if min_pep_len <= pep_len <= max_pep_len:
                    out_f.write(line)
