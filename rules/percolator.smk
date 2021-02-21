rule extract_pin:
    '''
        Extract compressed split PIN file.
    '''
    input:
        pin="out/{db}/split_pins/{ds}/{subset}/{sample}.{sdb}.pin.gz"
    output:
        pin=temp("out/{db}/split_pins/{ds}/{subset}/{sample}.{sdb}.pin")
    log:
        o="log/{db}/extract_pins/{ds}/{subset}/{sample}.{sdb}.out",
        e="log/{db}/extract_pins/{ds}/{subset}/{sample}.{sdb}.err"
    resources:
        mem = 8000
    threads: 1
    script:
        "../scripts/extract_pin.py"

rule percolator:
    '''
        Run the crux percolator algorithm to separate target from decoy matches.
    '''
    input:
        files=percolator_input_files
    output:
        file=temp("out/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}/percolator.target.psms.txt"),
        par="out/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}/percolator.params.txt",
        log="out/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}/percolator.log.txt"
    params:
        enzyme=percolator_enzyme
    log:
        o="log/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.out",
        e="log/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.err"
    conda:
        "../envs/crux.yaml"
    shadow: "minimal"
    threads: 1
    resources:
        mem = 34000,
        time = lambda wildcards, attempt: attempt * 480
    benchmark:
        "log/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.bmk"
    script:
        "../scripts/percolator.py"

