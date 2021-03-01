rule comet:
    '''
        Run comet on the files.
    '''
    input:
        data=comet_input_file,
        db="out/{db}/db/{td}.fasta",
        params=lambda wc: config["software"]["comet"]["params"].format(ds=wc.ds,subset=wc.subset)
    output:
        par="out/{db}/comet/{ds}/{subset}/{sample}.{td}/comet.params.txt",
        pin=temp("out/{db}/comet/{ds}/{subset}/{sample}.{td}/comet.target.pin"),
        log="out/{db}/comet/{ds}/{subset}/{sample}.{td}/comet.log.txt"
    log:
        o="log/{db}/comet/{ds}/{subset}/{sample}.{td}.out",
        e="log/{db}/comet/{ds}/{subset}/{sample}.{td}.err"
    benchmark:
        "log/{db}/comet/{ds}/{subset}/{sample}.{td}.bmk",
    threads: config["software"]["comet"]["threads"]
    conda:
        "../envs/crux.yaml"
    resources:
        mem = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 120
    script:
        "../scripts/comet.py"

rule split_pins:
    '''
        Split the comet output files by db. Merge targets and decoys.
    '''
    input:
        pins=["out/{{db}}/comet/{{ds}}/{{subset}}/{{sample}}.{td}/comet.target.pin".format(td=td)
              for td in ['target','decoy']]
    output:
        pin="out/{db}/split_pins/{ds}/{subset}/{sample}.{sdb}.pin.gz"
    log:
        o="log/{db}/split_pins/{ds}/{subset}/{sample}.{sdb}.out"
    resources:
        mem = 8000
    threads: 1
    run:
        import gzip
        with gzip.open(output.pin,"wt") as ofh:
            with open(input.pins[0]) as ifh: # write header from original pin
                ofh.write(ifh.readline())
            for f in input.pins:
                with open(f) as ifh:
                    for r in ifh:
                        if "#{}#".format(wildcards.sdb) in r:
                            ofh.write(r)
