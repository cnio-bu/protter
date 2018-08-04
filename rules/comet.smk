rule comet:
    '''
        Run comet on the files.
    '''
    input:
        data="res/data/prot/{ds}/{sample}.mzML",
        db="out/{db}/db/{td}.fasta"
    output:
        xml="out/{db}/comet/{ds}/{sample}.{td}.pep.xml"
    log:
        o="log/{db}/comet/{ds}/{sample}.{td}.out",
        e="log/{db}/comet/{ds}/{sample}.{td}.err"
    params:
        bin=config["software"]["search"]["comet"]["bin"],
        params=lambda wc: config["software"]["search"]["comet"]["params"].format(ds=wc.ds),
        basename=lambda wc: "out/{db}/comet/{ds}/{sample}.{td}".format(db=wc.db,sample=wc.sample,ds=wc.ds,td=wc.td)
    threads: config["software"]["search"]["comet"]["threads"]
    conda:
        "../envs/protter.yaml"
    resources:
        mem = 8000
    shell:"""
        export NSLOTS={threads}
        {params.bin} -P{params.params} -D{input.db} -N{params.basename} {input.data} > {log.o} 2> {log.e}
    """

