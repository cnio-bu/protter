rule msconvert:
    input:
        raw_file=msconvert_input_file
    output:
        mzml_file=msconvert_output_pattern()
    log:
        o="log/msconvert/{ds}/{sample}.out",
        e="log/msconvert/{ds}/{sample}.err"
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

        remote_host={config[software][msconvert][remote][hostname]}
        remote_user={config[software][msconvert][local][username]}
        remote_key={config[software][msconvert][local][private_key]}
        remote_addr="$remote_user@$remote_host"

        local_host={config[software][msconvert][local][hostname]}
        local_user={config[software][msconvert][local][username]}
        local_key={config[software][msconvert][local][private_key]}
        local_addr="$local_user@$local_host"

        pwiz_img={config[software][msconvert][docker_image]}
        raw_file="$(pwd)/{input.raw_file}"
        mzml_file="$(pwd)/{output.mzml_file}"

        ssh -i "$remote_key" "$remote_addr" 'bash -s' \
          < scripts/msconvert.sh -- -i "$local_key" -r "$local_addr" \
          "$pwiz_img" "$raw_file" "$mzml_file" 1>>{log.o} 2>>{log.e}

        echo [$(date --iso-8601=s)] {rule} done | tee -a {log.o} {log.e}
    """
