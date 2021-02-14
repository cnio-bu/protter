rule raw_to_mzml:
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
        ThermoRawFileParser -f=2 -m=0 --gzip \
            --input={input.raw_file} 1>{log.o} 2>{log.e}
    """
