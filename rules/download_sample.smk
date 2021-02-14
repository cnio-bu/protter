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
        wget --limit-rate="{resources.bandwidth}m" --tries 59 \
            --random-wait --wait=599 {params.file_url} \
            -O {output.data_file} 1>{log.o} 2>{log.e}

        rhash --printf="%h" {output.data_file} -o {output.cksum_file} \
            1>>{log.o} 2>>{log.e}

        exp_cksum={params.file_cksum}
        obs_cksum=$(cat {output.cksum_file})
        if [[ "$obs_cksum" != "$exp_cksum" ]]
        then
          echo "ERROR: checksum mismatch for file: {params.file_url}" >> {log.e}
          exit 1
        fi
    """
