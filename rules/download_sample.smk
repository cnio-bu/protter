rule download_sample:
    output:
        data_file=download_sample_output_pattern(),
        cksum_file=os.path.join(config["dataset_path"],"{ds}","{dl_file}.sha1")
    params:
        file_url=download_sample_file_url,
        file_cksum=download_sample_file_checksum
    log:
        o="log/download_sample/{ds}/{dl_file}.out",
        e="log/download_sample/{ds}/{dl_file}.err"
    benchmark:
        "log/download_sample/{ds}/{dl_file}.bmk",
    threads: 1
    conda:
        "../envs/environment.yaml"
    resources:
        mem = 1000,
        bandwidth = 2  # MB/s
    wildcard_constraints:
        gzip_ext = "(\.gz)?",
        fmt = "[^.]+"
    shell:"""
        bash scripts/download_file.sh -b {resources.bandwidth} \
          {params.file_url} {output.data_file} \
          1>{log.o} 2>{log.e}

        rhash --printf="%h" {output.data_file} -o {output.cksum_file}

        exp_cksum={params.file_cksum}
        obs_cksum=$(cat {output.cksum_file})
        if [[ "$obs_cksum" != "$exp_cksum" ]]
        then
          echo "ERROR: checksum mismatch for file: {params.file_url}"
          exit 1
        fi
    """
