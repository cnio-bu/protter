rule raw_to_mzml:
    '''
        Convert RAW file to mzML format.
    '''
    input:
        raw_file=raw_input_file
    output:
        mzml_file=mzml_output_pattern()
    log:
        o="log/raw_to_mzml/{ds}/{sample}.out",
        e="log/raw_to_mzml/{ds}/{sample}.err"
    benchmark:
        "log/raw_to_mzml/{ds}/{sample}.bmk",
    threads: 1
    conda:
        "../envs/raw_to_mzml.yaml"
    resources:
        mem = 8000
    shell:"""
        # We need to ensure logging is activated to make sure that
        # there is appropriate error output if anything goes wrong.
        ThermoRawFileParser -f=2 --gzip logging=1 \
            --input={input.raw_file} 1>{log.o} 2>{log.e}

        tmp_mzml_file=$(echo "{input.raw_file}" | sed 's/\.raw$//i')".mzML.gz"
        if [[ "$tmp_mzml_file" != "{output.mzml_file}" ]]
        then mv "$tmp_mzml_file" "{output.mzml_file}"
        fi
    """
