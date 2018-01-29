import os,glob

configfile: "config.json"

rule all:
    input:
        ["out/comet/{xp}".format(xp=os.path.splitext(os.path.basename(fname))[0]) for fname in glob.glob("res/data/raw/*.mzML")]

rule index_ML:
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

rule comet:
    input:
        conf=config["software"]["comet"]["conf"],
        raw="out/index_ML/{xp}.mzML",
        seq=config["seq_db"]
    output:
        basename="out/comet/{xp}"
    log:
        o="log/comet/{xp}.out",
        e="log/comet/{xp}.err"
    params:
        bin="bin/comet/{}".format(config["software"]["comet"]["bin"])
    threads: 32
    resources:
        mem = 8000
    shell:"""
        export NSLOTS={threads}
        {params.bin} -P{input.conf} -D{input.seq} -N{output.basename} {input.raw} > {log.o} 2> {log.e}
    """
