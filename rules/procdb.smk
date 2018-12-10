rule process_db:
    input:
        dbs=lambda wc: config["dbs"][wc.db]["paths"]
    output:
        db="out/{db}/db/target.fasta"
    conda:
        "../envs/protter.yaml"
    script:
        "scripts/process_db.py"

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
    params:
        tmp="out/{db}/db/decoyPYrat.tmp.fasta"
    conda:
        "../envs/protter.yaml"
    shell:"""
        python bin/decoyPYrat.py {input.db} -t {params.tmp}  -o {output.dec} --decoy_prefix decoy -k > {log.o} 2> {log.e}
        cat {input.db} {output.dec} >> {output.dbd} 2>> {log.e}
    """

