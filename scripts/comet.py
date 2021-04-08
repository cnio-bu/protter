from contextlib import redirect_stdout,redirect_stderr
import os
from subprocess import PIPE,run


def load_max_precursor_charge(comet_param_file):
    max_precursor_charge = 6  # default in Comet as of 2021-03-01

    with open(comet_param_file) as f:
        lines = f.readlines()

    # Comet accepts the value of the last occurrence of a parameter, so
    # we search the lines in reverse and stop as soon as we find it.
    for line in reversed(lines):
        content, *comments = line.split(sep="#",maxsplit=1)
        try:
            sep_idx = content.index("=")
        except ValueError:
            continue
        key = content[:sep_idx].strip()
        value = content[sep_idx + 1:].strip()
        if key == "max_precursor_charge":
            max_precursor_charge = int(value)
            break

    return max_precursor_charge


def write_percolator_header(max_precursor_charge,pin_file):
    '''
    Function to write a Percolator PIN file header.

    Adapted from the CometWritePercolator::WritePercolatorHeader
    method in the Comet MS/MS database search software under the
    following license. See also: comet-ms.sourceforge.net

        Copyright 2012 University of Washington

        Licensed under the Apache License, Version 2.0 (the "License");
        you may not use this file except in compliance with the License.
        You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

        Unless required by applicable law or agreed to in writing, software
        distributed under the License is distributed on an "AS IS" BASIS,
        WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
        See the License for the specific language governing permissions and
        limitations under the License.
    '''
    col_names = [
        "SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", "lnrSp", "deltLCn",
        "deltCn", "lnExpect", "Xcorr", "Sp", "IonFrac", "Mass", "PepLen"
    ] + [
        "Charge{}".format(i) for i in range(1, max_precursor_charge + 1)
    ] + [
        "enzN", "enzC", "enzInt", "lnNumSP", "dM", "absdM", "Peptide", "Proteins"
    ]

    with open(pin_file,"w") as f:
        header = "\t".join(col_names)
        f.write("{}\n".format(header))


data_file = snakemake.input.data
seqdb_file = snakemake.input.db
in_param_file = snakemake.input.param_file
out_param_file = snakemake.output.par
pin_file = snakemake.output.pin
comet_log_file = snakemake.output.log
out_log_file = snakemake.log.o
err_log_file = snakemake.log.e
threads = snakemake.threads


out_dir = os.path.dirname(pin_file)

with open(out_log_file,"w") as out_log_f, open(err_log_file,"w") as err_log_f:
    with redirect_stdout(out_log_f), redirect_stderr(err_log_f):
        cmd_args = [
            "crux","comet",
            "--num_threads",str(threads),
            "--verbosity","30",
            "--output_pepxmlfile","0",
            "--output_txtfile","0",
            "--output_percolatorfile","1",
            "--parameter-file",in_param_file,
            "--output-dir",out_dir,
            "--overwrite","T",
            data_file,
            seqdb_file
        ]
        run(cmd_args,stdout=PIPE,stderr=PIPE,check=True)

        # If Comet did not output a PIN file, check if this is due to a lack of searched
        # spectra in a job that was otherwise without error. If so, create a placeholder
        # PIN file to signal to Snakemake that the job was completed without error. Note
        # that this check currently depends on the Comet log verbosity being at least 30.
        if not os.path.isfile(pin_file):

            spectra_warning_found = False
            with open(comet_log_file) as f:
                for line in f:
                    line = line.rstrip()
                    if line.endswith("Warning - no spectra searched."):
                        spectra_warning_found = True

            if line.endswith("Return Code:0") and spectra_warning_found:
                max_precursor_charge = load_max_precursor_charge(in_param_file)
                write_percolator_header(max_precursor_charge, pin_file)
