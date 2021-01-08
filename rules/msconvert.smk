rule msconvert:
    input:
        raw_file=msconvert_input_file
    output:
        mzml_file=msconvert_output_pattern()
    params:
        rule_config=msconvert_rule_config()
    log:
        o="log/msconvert/{ds}/{sample}.out",
        e="log/msconvert/{ds}/{sample}.err"
    benchmark:
        "log/msconvert/{ds}/{sample}.bmk",
    threads: 1
    resources:
        mem = 8000,
        msconvert_slots = 1
    shell:"""
        echo [$(date --iso-8601=s)] {rule} starting | tee {log.o} {log.e}

        if [[ -z "{input.raw_file}" ]]
        then
          echo "cannot perform conversion - raw file not specified" >&2
          exit 1
        fi

        raw_file="$(pwd)/{input.raw_file}"
        mzml_file="$(pwd)/{output.mzml_file}"

        pwiz_image={params.rule_config[docker_image]}
        docker_host={params.rule_config[docker_host]}

        if [[ "$docker_host" == "local" ]]
        then
          bash scripts/msconvert.sh "$pwiz_image" \
            "$raw_file" "$mzml_file" 1>>{log.o} 2>>{log.e}
        else
          remote_addr={params.rule_config[remote_addr]}
          remote_key={params.rule_config[remote_key]}
          local_addr={params.rule_config[local_addr]}
          local_key={params.rule_config[local_key]}

          ssh -i "$local_key" "$remote_addr" 'bash -s' \
            < scripts/msconvert.sh -- -i "$remote_key" -r "$local_addr" \
            "$pwiz_image" "$raw_file" "$mzml_file" 1>>{log.o} 2>>{log.e}
        fi

        echo [$(date --iso-8601=s)] {rule} done | tee -a {log.o} {log.e}
    """
