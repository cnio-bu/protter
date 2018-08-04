rule percolator:
    '''
        Run the crux percolator agorithm to separate target from decoy matches.
    '''
    input:
        targets = lambda wc: ["out/{db}/{sw}/{ds}/{sample}.target.pep.xml".format(db=wc.db,sw=wc.sw,ds=wc.ds,sample=sample) for sample in get_samples(wc.ds,wc.grouping)[wc.group]],
        decoys = lambda wc: ["out/{db}/{sw}/{ds}/{sample}.decoy.pep.xml".format(db=wc.db,sw=wc.sw,ds=wc.ds,sample=sample) for sample in get_samples(wc.ds,wc.grouping)[wc.group]],
    output:
        d="out/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}",
        f="out/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}/percolator.target.peptides.txt",
        l="out/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}/input_list",
        p=temp("out/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}/make-pin.pin")
    log:
        o="log/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}.out",
        e="log/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}.err"
    threads: 1
    conda:
        "../envs/protter.yaml"
    resources:
        mem = 34000
    benchmark:
        "log/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}.bmk"
    run:
        if wildcards.em == "tdc":
            tdc = "T"
        elif wildcards.em == "mix-max":
            tdc = "F"
        else:
            raise ValueError("Wrong value for --tdc: {}".format(wildcards.em))
        with open(output.l,'w') as ofh:
            for t in input.targets:
                ofh.write("{}\n".format(t))
            for t in input.decoys:
                ofh.write("{}\n".format(t))
        shell("bin/crux/crux percolator --tdc {tdc} --overwrite T --protein T --fido-empirical-protein-q T --output-dir {output.d} --list-of-files T {output.l} > {log.o} 2> {log.e}")

