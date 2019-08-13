import subprocess,shlex

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

with open(log.o,'w') as lfh:
    with open(log.e,'w') as efh:
        cmd = "crux percolator --decoy_prefix decoy_ --tdc {tdc} --overwrite T --protein T --fido-empirical-protein-q T --output-dir {d} {i}".format(tdc = tdc,d=par.d,i=inp.files)
        subprocess.call(shlex.split(cmd),stdout=lfh,stderr=efh)
