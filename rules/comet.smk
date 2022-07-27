rule comet:
    '''
        Run comet on the files.
    '''
    input:
        data=comet_input_file,
        db="out/{db}/db/{td}.fasta",
        param_file=lambda wc: config["software"]["comet"]["param_file"].format(ds=wc.ds,
                                                                               subset=wc.subset)
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
        mem_mb = lambda wildcards, attempt: attempt * 12000,
        time = lambda wildcards, attempt: attempt * 120
    script:
        "../scripts/comet.py"

rule filter_pins:
    '''
        Filter comet output files by peptide length.

        Comet release 2020.01 supports a 'peptide_length_range' parameter.
        As of March 2021, the version of Crux toolkit available on Bioconda
        (Crux version 3.2-0d57cff) includes Comet version 2018.01 rev. 0,
        which does not support this parameter. As a stopgap, this rule
        filters PIN files by peptide length.
    '''
    input:
        pin="out/{db}/comet/{ds}/{subset}/{sample}.{td}/comet.target.pin"
    output:
        pin=temp("out/{db}/filter_pins/{ds}/{subset}/{sample}.{td}/comet.target.pin")
    params:
        peptide_length_range=config["software"]["comet"]["params"]["peptide_length_range"]
    log:
        o="log/{db}/filter_pins/{ds}/{subset}/{sample}.{td}.out",
        e="log/{db}/filter_pins/{ds}/{subset}/{sample}.{td}.err"
    threads: 1
    resources:
        mem_mb = 1000
    script:
        "../scripts/filter_pins.py"

rule split_pins:
    '''
        Split the comet output files by db. Merge targets and decoys.
    '''
    input:
        pins=["out/{{db}}/filter_pins/{{ds}}/{{subset}}/{{sample}}.{td}/comet.target.pin".format(td=td)
              for td in ['target','decoy']]
    output:
        pin="out/{db}/split_pins/{ds}/{subset}/{sample}.{sdb}.pin.gz"
    log:
        o="log/{db}/split_pins/{ds}/{subset}/{sample}.{sdb}.out"
    resources:
        mem_mb = 8000
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
