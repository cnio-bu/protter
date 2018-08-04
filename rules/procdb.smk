rule process_db:
    input:
        dbs=lambda wc: config["dbs"][wc.db]["paths"]
    output:
        db="out/{db}/db/target.fasta"
    conda:
        "../envs/protter.yaml"
    run:
        from tstk.io import parsefastx

        with open(output.db,"w") as ofh:
            for fpath in input.dbs:
                with open(fpath) as fh:
                    for r in parsefastx(fh):
                        ofh.write(">{}\n{}\n".format(r[0],r[1].replace("L","I")))

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

