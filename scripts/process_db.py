from tstk.io import parsefastx

with open(snakemake.output.db,"w") as ofh:
    for fpath in snakemake.input.dbs:
        with open(fpath) as fh:
            for r in parsefastx(fh):
                ofh.write(">{}\n{}\n".format(r[0],r[1].replace("L","I")))
