rule group_pins:
    '''
        Merge PIN files by group.
    '''
    # Though it can accept multiple input files, Percolator has raised filestream errors
    # when run as a shadow rule with multiple PIN files, so we merge them here first.
    input:
        pins=group_pin_files
    output:
        pin=temp("out/{db}/group_pins/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}.pin")
    log:
        o="log/{db}/group_pins/{ds}/{subset}/{em}/{grouping}.{group}.{sdb}.out",
        e="log/{db}/group_pins/{ds}/{subset}/{em}/{grouping}.{group}.{sdb}.err"
    resources:
        mem = 8000
    threads: 1
    script:
        "../scripts/group_pins.py"

rule percolator:
    '''
        Run the crux percolator algorithm to separate target from decoy matches.
    '''
    input:
        file="out/{db}/group_pins/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}.pin"
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
        mem = lambda wildcards, attempt: attempt * 64000,
        time = 1440
    benchmark:
        "log/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.bmk"
    script:
        "../scripts/percolator.py"
