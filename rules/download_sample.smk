rule download_sample:
    output:
        data_file=download_sample_output_pattern(),
        cksum_file=os.path.join(config["download_path"],"{ds}","{dl_file}.sha1")
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
        "../envs/download.yaml"
    resources:
        mem = 1000,
        bandwidth = 2  # MB/s
    shell:"""
        curl --silent --show-error --retry 17 --limit-rate "{resources.bandwidth}m" \
            {params.file_url} --output {output.data_file} 1>{log.o} 2>{log.e}

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


rule confirm_checksums:
    input:
        cksum_files=lambda wc: [
            os.path.join(config["download_path"],wc.ds,"{}.sha1".format(dl_file))
            for dl_file in downloads.loc[pd.IndexSlice[wc.ds],"dl_file"]
        ],
        config_file="config.yaml"
    output:
        cksum_file=os.path.join(config["download_path"],"{ds}","sha1checksums.txt")
    log:
        o="log/download_sample/{ds}/confirm_checksums.out",
        e="log/download_sample/{ds}/confirm_checksums.err"
    threads: 1
    conda:
        "../envs/common.yaml"
    script:
        "../scripts/confirm_checksums.py"
