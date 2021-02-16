rule percolator:
    '''
        Run the crux percolator agorithm to separate target from decoy matches.
    '''
    input:
        files=percolator_input_files
    output:
        f="out/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}/percolator.target.psms.txt"
    params:
        d="out/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}",
        enzyme=percolator_enzyme
    log:
        o="log/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.out",
        e="log/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.err"
    conda:
        "../envs/crux.yaml"
    threads: 1
    resources:
        mem = 34000,
        time = lambda wildcards, attempt: attempt * 480
    benchmark:
        "log/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.bmk"
    script:
        "../scripts/percolator.py"

