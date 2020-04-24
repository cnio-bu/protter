
import json


def load_sample_metadata(wildcards,input):
    # Wildcards variable is not used, but is required by Snakemake.
    with open(input.proxy_file) as f:
        return json.load(f)


rule download_sample:
    input:
        proxy_file=config["dataset_path"] + "/{ds}/{sample}.{fmt}{gzip_ext}_proxy.json"
    output:
        data_file=download_sample_output_pattern()
    params:
        meta=load_sample_metadata
    log:
        o="log/download_sample/{ds}/{sample}.{fmt}{gzip_ext}.out",
        e="log/download_sample/{ds}/{sample}.{fmt}{gzip_ext}.err"
    benchmark:
        "log/download_sample/{ds}/{sample}.{fmt}{gzip_ext}.bmk",
    threads: 1
    resources:
        mem = 1000,
        bandwidth = 2  # MB/s
    wildcard_constraints:
        gzip_ext = "(\.gz)?"
    shell:"""
        bash scripts/download_file.sh -b {resources.bandwidth} \
          {params[meta][url]} {output.data_file} \
          1>{log.o} 2>{log.e}
    """
