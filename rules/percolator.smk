rule percolator:
    '''
        Run the crux percolator agorithm to separate target from decoy matches.
    '''
    input:
        files = lambda wc: ["out/{db}/split_pins/{ds}/{sample}.{sdb}.pin".format(db=wc.db,ds=wc.ds,sample=sample,sdb=wc.sdb) for sample in get_samples(wc.ds,wc.grouping)[wc.group]]
    output:
        f="out/{db}/percolator/{ds}/{em}/{grouping}/{group}/{sdb}/percolator.target.peptides.txt"
    params:
        d="out/{db}/percolator/{ds}/{em}/{grouping}/{group}/{sdb}"
    log:
        o="log/{db}/percolator/{ds}/{em}/{grouping}/{group}.{sdb}.out",
        e="log/{db}/percolator/{ds}/{em}/{grouping}/{group}.{sdb}.err"
    conda:
        "../envs/environment.yaml"
    threads: 1
    resources:
        mem = 34000
    benchmark:
        "log/{db}/percolator/{ds}/{em}/{grouping}/{group}.{sdb}.bmk"
    script:
        "../scripts/percolator.py"

