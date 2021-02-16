rule process_psms:
    '''
        Post-process output PSM data.
    '''
    input:
        psms="out/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}/percolator.target.psms.txt",
        target_meta_file="out/{db}/db/target_metadata.tsv",
        sample_meta_file=config["samples"]
    output:
        psms="out/{db}/percolator/{ds}/{subset}/{em}/{grouping}/{group}/{sdb}/percolator.target.psms.proc.txt"
    params:
        pep_cutoff=config["software"]["percolator"]["pep_cutoff"]
    log:
        o="log/{db}/process_psms/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.out",
        e="log/{db}/process_psms/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.err"
    benchmark:
        "log/{db}/process_psms/{ds}/{subset}/{em}/{grouping}/{group}.{sdb}.bmk"
    threads: 1
    conda:
        "../envs/common.yaml"
    resources:
        mem = 8000
    script:
        "../scripts/process_psms.py"