import os,glob

configfile: "config.json"

rule all:
    input:
        ["out/comet/{xp}".format(xp=os.path.splitext(os.path.basename(fname))[0]) for fname in glob.glob("res/data/raw/*.mzML")]

rule ML2XML:
    input:
        ML="res/data/raw/{xp}.mzML",
    output:
        XML="out/ML2XML/{xp}.mzXML"
    log:
        o="log/ML2XML/{xp}.out",
        e="log/ML2XML/{xp}.err"
    threads: 1
    resources:
        mem = 1000
    shell:"""
        msconvert {input.ML} --mzXML -o out/ML2XML > {log.o} 2> {log.e}
    """

rule comet:
    input:
        conf=config["software"]["comet"]["conf"],
        raw="out/ML2XML/{xp}.mzXML",
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
