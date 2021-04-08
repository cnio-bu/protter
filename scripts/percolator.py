from contextlib import redirect_stdout,redirect_stderr
from itertools import islice
import os
from subprocess import PIPE,run


def create_empty_file(file_path):
    with open(file_path,"w") as f:
        pass


def write_default_psm_file(psm_file):
    """Write a default Percolator PSM file."""

    col_names = [
        "file_idx", "scan", "charge", "spectrum precursor m/z",
        "spectrum neutral mass", "peptide mass", "percolator score",
        "percolator q-value", "percolator PEP", "total matches/spectrum",
        "sequence", "protein id", "flanking aa"
    ]

    with open(psm_file,"w") as f:
        header = "\t".join(col_names)
        f.write("{}\n".format(header))


wc = snakemake.wildcards
par = snakemake.params
out = snakemake.output
inp = snakemake.input
log = snakemake.log

if wc.em == "tdc":
    tdc = "T"
elif wc.em == "mix-max":
    tdc = "F"
else:
    raise ValueError("Wrong value for --tdc: {}".format(wc.em))

out_dir = os.path.dirname(out.file)

subset_max_train = 1000000

with open(log.o,"w") as out_log_f, open(log.e,"w") as err_log_f:
    with redirect_stdout(out_log_f), redirect_stderr(err_log_f):

        # Short-circuit for PIN file lacking data.
        with open(inp.file) as f:
            header_line, *data_lines = islice(f,2)
            if not data_lines:
                write_default_psm_file(out.file)
                create_empty_file(out.par)
                create_empty_file(out.log)
                sys.exit(0)

        cmd_args = [
            "crux","percolator",
            "--decoy_prefix","decoy_",
            "--tdc",tdc,
            "--overwrite","T",
            "--protein","T",
            "--fido-empirical-protein-q","T",
            "--subset-max-train",str(subset_max_train),
            "--protein-enzyme",par.enzyme,
            "--output-dir",out_dir,
            inp.file
        ]
        run(cmd_args,stdout=PIPE,stderr=PIPE,check=True)
