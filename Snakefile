import os,glob
from tstk.io import parsefastx
import xml.etree.ElementTree as ET

configfile: "config.yaml"

def get_samples(ds):
    fpath = "res/data/prot/{ds}".format(ds=ds)
    csvpath  = fpath + "/groups.csv"
    groups = {}
    with open(csvpath) as fh:
        for r in fh:
            r = r.rstrip("\n")
            g,f = r.split(",")
            f = f.rstrip(".mzML")
            try:
                groups[g].append(f)
            except KeyError:
                groups[g] = [f]
    return groups

def input_all():
    for software in config["software"]["search"]:
        if config["software"]["search"][software]["enabled"]:
            for ds in config["datasets"]:
                if config["datasets"][ds]["enabled"]:
                    for group in get_samples(ds):
                        for db in config["dbs"]:
                            if config["dbs"][db]["enabled"]:
                                yield "out/{db}/percolator/{sw}/{ds}/{group}/percolator.target.peptides.txt".format(group=group,sw=software,db=db,ds=ds)

rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        input_all()

rule convert_leucines:
    input:
        db=lambda wc: config["dbs"][wc.db]["path"]
    output:
        db="out/{db}/db/target.fasta"
    run:
        with open(input.db) as fh:
            with open(output.db,"w") as ofh:
                for r in parsefastx(fh):
                    ofh.write(">{}\n{}\n".format(r[0],r[1].replace("L","I")))

#rule raw2mzML:
#    input:
#        raw="res/data/prot/{ds}/{sample}.raw"
#    output:
#        mzML="res/data/prot/{ds}/{sample}.mzML"
#    log:    
#        o="log/raw2mzML/{ds}/{sample}.out",
#        e="log/raw2mzML/{ds}/{sample}.err"
#    shell:'''
#        msconvert {input.raw} --mzXML --outfile {output.mzML} > {log.o} 2> {log.e}
#    '''
#
#rule raw2mgf:
#    input:
#        raw="res/data/prot/{ds}/{sample}.raw"
#    output:
#        mgf="res/data/prot/{ds}/{sample}.mgf"
#    log:    
#        o="log/raw2mgf/{ds}/{sample}.out",
#        e="log/raw2mgf/{ds}/{sample}.err"
#    shell:'''
#        msconvert {input.raw} --mgf --outfile {output.mgf} > {log.o} 2> {log.e}
#    '''

rule add_decoys:
    '''
        Use decoPYrat to add decoy sequences to the original database.
    '''
    input:
        db="out/{db}/db/target.fasta"
    output:
        dec="out/{db}/db/decoy.fasta",
        dbd="out/{db}/db/target_and_decoy.fasta",
    log:
        o="log/add_decoys/{db}.out",
        e="log/add_decoys/{db}.err"
    params:
        tmp="out/{db}/db/decoyPYrat.tmp.fasta"
    shell:"""
        python bin/decoyPYrat.py {input.db} -t {params.tmp}  -o {output.dec} --decoy_prefix decoy > {log.o} 2> {log.e}
        cat {input.db} {output.dec} >> {output.dbd} 2>> {log.e}
    """

rule comet:
    '''
        Run comet on the files against the "decoyed" database.
    '''
    input:
        data="res/data/prot/{ds}/{sample}.mzML",
        db="out/{db}/db/{td}.fasta",
    output:
        xml="out/{db}/comet/{ds}/{sample}.{td}.pep.xml"
    log:
        o="log/{db}/comet/{ds}/{sample}.{td}.out",
        e="log/{db}/comet/{ds}/{sample}.{td}.err"
    params:
        bin=config["software"]["search"]["comet"]["bin"],
        params=lambda wc: config["software"]["search"]["comet"]["params"].format(ds=wc.ds),
        basename=lambda wc: "out/{db}/comet/{ds}/{sample}.{td}".format(db=wc.db,sample=wc.sample,ds=wc.ds,td=wc.td)
    threads: config["software"]["search"]["comet"]["threads"]
    resources:
        mem = 8000
    shell:"""
        export NSLOTS={threads}
        {params.bin} -P{params.params} -D{input.db} -N{params.basename} {input.data} > {log.o} 2> {log.e}
    """

rule xtandem_taxonomy:
    input:
        db="out/{db}/db/{td}.fasta"
    output:
        tax="out/{db}/xtandem/taxonomy/{td}.xtandem.taxonomy.xml"
    run:
        root = ET.Element("biomla", label="x! taxon-to-file matching list")
        taxon = ET.SubElement(root,"taxon", label=config["dbs"][wildcards.db]["tax"])
        f = ET.SubElement(taxon,"file", format="peptide", URL=input.db)

        tree = ET.ElementTree(root)
        tree.write(output.tax, xml_declaration=True)

rule xtandem:
    '''
        Run xtandem on the files against the "decoyed" database.
    '''
    input:
        data="res/data/prot/{ds}/{sample}.mzML",
        tax="out/{db}/xtandem/taxonomy/{td}.xtandem.taxonomy.xml"
    output:
        xml=touch("out/{db}/xtandem/{ds}/{sample}.{td}.t.xml"),
        conf="out/{db}/xtandem/{ds}_conf/{sample}.{td}.xml"
    log:
        o="log/{db}/xtandem/{ds}/{sample}.{td}.out",
        e="log/{db}/xtandem/{ds}/{sample}.{td}.err"
    params:
        bin=config["software"]["search"]["xtandem"]["bin"],
        params=lambda wc: config["software"]["search"]["xtandem"]["params"].format(ds=wc.ds),
        basename=lambda wc: "out/{db}/xtandem/{ds}/{sample}.{td}".format(db=wc.db,ds=wc.ds,sample=wc.sample,td=wc.td)
    threads: config["software"]["search"]["xtandem"]["threads"]
    resources:
        mem = 8000
    run:
        shell("rm -f {params.basename}*.xml")

        tree = ET.parse("res/conf/xtandem/xtandem.input_template.xml")
        root = tree.getroot()

        for n in root.iter("note"):
            try:
                l = n.attrib["label"]
            except:
                continue
            
            if l == "list path, default parameters":
                n.text = params.params
            elif l == "list path, taxonomy information":
                n.text = input.tax
            elif l == "protein, taxon":
                n.text = config["dbs"][wildcards.db]["tax"]
            elif l == "spectrum, path":
                n.text = input.data
            elif l == "output, path":
                n.text = output.xml

        tree.write(output.conf, xml_declaration=True)

        shell("{params.bin} {output.conf} > {log.o} 2> {log.e}")
        shell("mv {params.basename}*.xml {output.xml}")

rule xtandem2pepxml:
    '''
        Convert t.xml files to pep.xml
    '''
    input:
        txml="out/{db}/xtandem/{ds}/{sample}.{td}.t.xml",
    output:
        pepxml="out/{db}/xtandem/{ds}/{sample}.{td}.pep.xml",
    log:
        o="log/{db}/xtandem2pepxml/{ds}/{sample}.{td}.out",
        e="log/{db}/xtandem2pepxml/{ds}/{sample}.{td}.err"
    params:
    threads: 1
    resources:
        mem = 4000
    shell:'''
        pepxmltk.py {input.txml} {output.pepxml}   
    '''

rule percolator:
    '''
        Run the crux percolator agorithm to separate target from decoy matches.
    '''
    input:
        targets = lambda wc: ["out/{db}/{sw}/{ds}/{sample}.target.pep.xml".format(db=wc.db,sw=wc.sw,ds=wc.ds,sample=sample) for sample in get_samples(wc.ds)[wc.group]],
        decoys = lambda wc: ["out/{db}/{sw}/{ds}/{sample}.decoy.pep.xml".format(db=wc.db,sw=wc.sw,ds=wc.ds,sample=sample) for sample in get_samples(wc.ds)[wc.group]]
    output:
        d="out/{db}/percolator/{sw}/{ds}/{group}",
        f="out/{db}/percolator/{sw}/{ds}/{group}/percolator.target.peptides.txt",
        l="out/{db}/percolator/{sw}/{ds}/{group}/input_list"
    log:
        o="log/{db}/percolator/{sw}/{ds}/{group}.out",
        e="log/{db}/percolator/{sw}/{ds}/{group}.err"
    threads: 1
    resources:
        mem = 8000
    run:
        with open(output.l,'w') as ofh:
            for t in input.targets:
                ofh.write("{}\n".format(t))
        shell("bin/crux/crux percolator --overwrite T --protein T --fido-empirical-protein-q T --output-dir {output.d} --list-of-files T {output.l} > {log.o} 2> {log.e}")
