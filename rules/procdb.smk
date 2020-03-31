rule process_db:
    input:
        paths=lambda wc: [config["dbs"][wc.db]["paths"][k] for k in config["dbs"][wc.db]["paths"]]
    params:
        dbs=lambda wc: [k for k in config["dbs"][wc.db]["paths"]]
    output:
        db="out/{db}/db/target.fasta"
    conda:
        "../envs/environment.yaml"
    log:
        o="log/process_db/{db}.out"
    benchmark:
        "log/process_db/{db}.bmk",
    script:
        "../scripts/process_db.py"

rule add_decoys:
    '''
        Use decoPYrat to generate decoy sequences.
    '''
    input:
        db="out/{db}/db/target.fasta"
    output:
        dec="out/{db}/db/decoy.fasta",
        dbd="out/{db}/db/target_and_decoy.fasta",
    log:
        o="log/add_decoys/{db}.out",
        e="log/add_decoys/{db}.err"
    benchmark:
        "log/add_decoys/{db}.bmk",
    params:
        tmp="out/{db}/db/decoyPYrat.tmp.fasta"
    conda:
        "../envs/environment.yaml"
    shell:"""
        decoypyrat {input.db} -t {params.tmp}  -o {output.dec} --decoy_prefix decoy -k > {log.o} 2> {log.e}
        cat {input.db} {output.dec} >> {output.dbd} 2>> {log.e}
    """

