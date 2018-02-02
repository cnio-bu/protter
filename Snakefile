import os,glob
from tstk.io import parsefastx

configfile: "config.yaml"

rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        ["out/dbs/{db}/comet/{xp}.pep.xml".format(db=db,xp=os.path.splitext(os.path.basename(fname))[0]) for db in config["active_dbs"] for fname in glob.glob("res/data/raw/*.mzML")]

rule index_ML:
    '''
        Use msconvert to generate indexed versions of the original .mzML files.
    '''
    input:
        ML="res/data/raw/{xp}.mzML",
    output:
        ML="out/index_ML/{xp}.mzML"
    log:
        o="log/index_ML/{xp}.out",
        e="log/index_ML/{xp}.err"
    threads: 1
    resources:
        mem = 1000
    shell:"""
        msconvert {input.ML} --mzML -o out/index_ML > {log.o} 2> {log.e}
    """

rule convert_leucines:
    input:
        db=lambda wc: config["dbs"][wc.db]["path"]
    output:
        db="out/dbs/{db}/convert_leucines/db.fasta"
    run:
        with open(input.db) as fh:
            with open(output.db,"w") as ofh:
                for r in parsefastx(fh):
                    ofh.write(">{}\n{}\n".format(r[0],r[1].replace("L","I")))

rule add_decoys:
    '''
        Use decyPYrat to add decoy sequences to the original database.
    '''
    input:
        db="out/dbs/{db}/convert_leucines/db.fasta"
    output:
        dec="out/dbs/{db}/add_decoys/decoys.fasta",
        dbd="out/dbs/{db}/add_decoys/db_and_decoys.fasta",
    log:
        o="log/dbs/{db}/add_decoys/{db}.out",
        e="log/dbs/{db}/add_decoys/{db}.err"
    params:
        tmp="out/dbs/{db}/add_decoys/decoyPYrat.tmp.fasta"
    shell:"""
        python bin/decoyPYrat.py {input.db} -t {params.tmp}  -o {output.dec} --decoy_prefix decoy > {log.o} 2> {log.e}
        cat {input.db} {output.dec} >> {output.dbd} 2> {log.e}
    """

rule comet:
    '''
        Run comet on the files against the "decoyed" database.
    '''
    input:
        conf=config["software"]["comet"]["conf"],
        raw="out/index_ML/{xp}.mzML",
        seq="out/dbs/{db}/add_decoys/db_and_decoys.fasta"
    output:
        xml="out/dbs/{db}/comet/{xp}.pep.xml"
    log:
        o="log/dbs/{db}/comet/{xp}.out",
        e="log/dbs/{db}/comet/{xp}.err"
    params:
        bin="bin/comet/{}".format(config["software"]["comet"]["bin"]),
        basename=lambda wc: "out/comet/{xp}".format(xp=wc.xp)
    threads: config["software"]["comet"]["threads"]
    resources:
        mem = 8000
    shell:"""
        export NSLOTS={threads}
        {params.bin} -P{input.conf} -D{input.seq} -N{params.basename} {input.raw} > {log.o} 2> {log.e}
    """
