rule comet:
    '''
        Run comet on the files.
    '''
    input:
        data=comet_input_file,
        db="out/{db}/db/{td}.fasta",
        params=lambda wc: config["software"]["comet"]["params"].format(ds=wc.ds)
    output:
        par="out/{db}/comet/{ds}/{sample}.{td}/comet.params.txt",
        pin="out/{db}/comet/{ds}/{sample}.{td}/comet.target.pin",
        log="out/{db}/comet/{ds}/{sample}.{td}/comet.log.txt"
    log:
        o="log/{db}/comet/{ds}/{sample}.{td}.out",
        e="log/{db}/comet/{ds}/{sample}.{td}.err"
    benchmark:
        "benchmark/{db}/comet/{ds}/{sample}.{td}.tsv",
    params:
        basename="out/{db}/comet/{ds}/{sample}.{td}"
    threads: config["software"]["comet"]["threads"]
    conda:
        "../envs/environment.yaml"
    resources:
        mem = lambda wildcards, attempt: attempt * 8000
    shell:"""
        crux comet --num_threads {threads} --output_pepxmlfile 0 --output_txtfile 0 --output_percolatorfile 1 --parameter-file {input.params} --output-dir {params.basename} --overwrite T {input.data} {input.db} > {log.o} 2> {log.e}
    """

rule split_pins:
    '''
        Split the comet output files by db. Merge targets and decoys.
    '''
    input:
        pins=["out/{{db}}/comet/{{ds}}/{{sample}}.{td}/comet.target.pin".format(td=td) for td in ['target','decoy']]
    output:
        pin="out/{db}/split_pins/{ds}/{sample}.{sdb}.pin"
    log:
        o="log/{db}/split_pins/{ds}/{sample}.{sdb}.out"
    resources:
        mem = 8000
    threads: 1
    run:
        with open(output.pin,"w") as ofh:
            with open(input.pins[0]) as ifh: # write header from original pin
                ofh.write(ifh.readline())
            for f in input.pins:
                with open(f) as ifh:
                    for r in ifh:
                        if "#{}#".format(wildcards.sdb) in r:
                            ofh.write(r)
