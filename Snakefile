import os,glob
from tstk.io import parsefastx

configfile: "config.yaml"

rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        ["out/{db}/comet/{xp}.pep.xml".format(db=db,xp=os.path.splitext(os.path.basename(fname))[0]) for db in config["active_dbs"] for fname in glob.glob("res/data/raw/*.mzML")],
        ["out/common/raw/{xp}.mgf".format(xp=os.path.splitext(os.path.basename(fname))[0]) for fname in glob.glob("res/data/raw/*.mzML")]

rule index_mzML:
    '''
        Use msconvert to generate indexed versions of the original .mzML files.
    '''
    input:
        ML="res/data/raw/{xp}.mzML",
    output:
        ML="out/common/raw/{xp}.mzML"
    log:
        o="log/index_mzML/{xp}.out",
        e="log/index_mzML/{xp}.err"
    threads: 1
    resources:
        mem = 1000
    shell:"""
        msconvert {input.ML} --mzML -o out/common/raw > {log.o} 2> {log.e}
    """

rule convert_mgf:
    '''
        Use msconvert to convert mzML files to mgf (for X!Tandem).
    '''
    input:
        ML="res/data/raw/{xp}.mzML",
    output:
        ML="out/common/raw/{xp}.mgf"
    log:
        o="log/convert_mgf/{xp}.out",
        e="log/convert_mgf/{xp}.err"
    threads: 1
    resources:
        mem = 1000
    shell:"""
        msconvert {input.ML} --mgf -o out/common/raw > {log.o} 2> {log.e}
    """

rule convert_leucines:
    input:
        db=lambda wc: config["dbs"][wc.db]["path"]
    output:
        db="out/{db}/db/db.fasta"
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
        db="out/{db}/db/db.fasta"
    output:
        dec="out/{db}/db/decoys.fasta",
        dbd="out/{db}/db/db_and_decoys.fasta",
    log:
        o="log/add_decoys/{db}.out",
        e="log/add_decoys/{db}.err"
    params:
        tmp="out/{db}/db/decoyPYrat.tmp.fasta"
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
        raw="out/common/raw/{xp}.mzML",
        db="out/{db}/db/db_and_decoys.fasta",
    output:
        xml="out/{db}/comet/{xp}.pep.xml"
    log:
        o="log/{db}/comet/{xp}.out",
        e="log/{db}/comet/{xp}.err"
    params:
        bin="bin/comet/{}".format(config["software"]["comet"]["bin"]),
        basename=lambda wc: "out/{db}/comet/{xp}".format(db=wc.db,xp=wc.xp)
    threads: config["software"]["comet"]["threads"]
    resources:
        mem = 8000
    shell:"""
        export NSLOTS={threads}
        {params.bin} -P{input.conf} -D{input.db} -N{params.basename} {input.raw} > {log.o} 2> {log.e}
    """
