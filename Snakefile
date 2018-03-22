import os,glob
from tstk.io import parsefastx
import xml.etree.ElementTree as ET

configfile: "config.yaml"

def input_files(software,iftype,oftype):
    if config["software"][software]["enabled"]:
        for ds in config["datasets"]:
            if config["datasets"][ds]["enabled"]:
                for db in config["active_dbs"]:
                    for fname in glob.glob("res/data/prot/{iftype}/{ds}/*.{iftype}".format(ds=ds,iftype=iftype)):
                        #for sdb in ['db','decoys','db_and_decoys']:
                        for sdb in ['db','decoys']:
                            yield "out/{db}/{sw}/{ds}/{sdb}/{xp}.{oftype}".format(sdb=sdb,sw=software,db=db,xp=os.path.splitext(os.path.basename(fname))[0],ds=ds,oftype=oftype)
rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        input_files("comet","mzML","pep.xml"),
        input_files("comet","raw","pep.xml"),
        input_files("xtandem","mgf","pep.xml"),
        input_files("xtandem","raw","pep.xml")

rule convert_leucines:
    input:
        db=lambda wc: config["dbs"][wc.db]["path"]
    output:
        db="out/{db}/db/db.fasta"
    run:
        with open(input.db) as fh:
            with open(output.db,"w") as ofh:
                for r in parsefastx(fh):
                    ofh.write(">{}\n{}\n".format(r[0],r[1].replace("L","I")))

rule raw2mzML:
    input:
        raw="res/data/prot/raw/{ds}/{xp}.raw"
    output:
        mzML="res/data/prot/mzML/{ds}/{xp}.mzML"
    log:    
        o="log/raw2mzML/{ds}/{xp}.out",
        e="log/raw2mzML/{ds}/{xp}.err"
    shell:'''
        msconvert {input.raw} --mzXML --outfile {output.mzML} > {log.o} 2> {log.e}
    '''

rule raw2mgf:
    input:
        raw="res/data/prot/raw/{ds}/{xp}.raw"
    output:
        mgf="res/data/prot/mgf/{ds}/{xp}.mgf"
    log:    
        o="log/raw2mgf/{ds}/{xp}.out",
        e="log/raw2mgf/{ds}/{xp}.err"
    shell:'''
        msconvert {input.raw} --mgf --outfile {output.mgf} > {log.o} 2> {log.e}
    '''

rule add_decoys:
    '''
        Use decoPYrat to add decoy sequences to the original database.
    '''
    input:
        db="out/{db}/db/db.fasta"
    output:
        dec="out/{db}/db/decoys.fasta",
        dbd="out/{db}/db/db_and_decoys.fasta",
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
        data="res/data/prot/mzML/{ds}/{xp}.mzML",
        db="out/{db}/db/{sdb}.fasta",
    output:
        xml="out/{db}/comet/{ds}/{sdb}/{xp}.pep.xml"
    log:
        o="log/{db}/comet/{ds}/{sdb}/{xp}.out",
        e="log/{db}/comet/{ds}/{sdb}/{xp}.err"
    params:
        bin=config["software"]["comet"]["bin"],
        params=lambda wc: config["software"]["comet"]["params"].format(ds=wc.ds),
        basename=lambda wc: "out/{db}/comet/{ds}/{sdb}/{xp}".format(db=wc.db,xp=wc.xp,ds=wc.ds,sdb=wc.sdb)
    threads: config["software"]["comet"]["threads"]
    resources:
        mem = 8000
    shell:"""
        export NSLOTS={threads}
        {params.bin} -P{params.params} -D{input.db} -N{params.basename} {input.data} > {log.o} 2> {log.e}
    """

rule xtandem_taxonomy:
    input:
        db="out/{db}/db/{sdb}.fasta"
    output:
        tax="out/{db}/xtandem/taxonomy/{sdb}.xtandem.taxonomy.xml"
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
        data="res/data/prot/mgf/{ds}/{xp}.mgf",
        tax="out/{db}/xtandem/taxonomy/{sdb}.xtandem.taxonomy.xml"
    output:
        xml="out/{db}/xtandem/{ds}/{sdb}/{xp}.t.xml",
        conf="out/{db}/xtandem/{ds}_conf/{sdb}/{xp}.xml"
    log:
        o="log/{db}/xtandem/{ds}/{sdb}/{xp}.out",
        e="log/{db}/xtandem/{ds}/{sdb}/{xp}.err"
    params:
        bin=config["software"]["xtandem"]["bin"],
        params=lambda wc: config["software"]["xtandem"]["params"].format(ds=wc.ds),
        basename=lambda wc: "out/{db}/xtandem/{ds}/{sdb}/{xp}".format(db=wc.db,ds=wc.ds,xp=wc.xp,sdb=wc.sdb)
    threads: config["software"]["xtandem"]["threads"]
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
        txml="out/{db}/xtandem/{ds}/{sdb}/{xp}.t.xml",
    output:
        pepxml="out/{db}/xtandem/{ds}/{sdb}/{xp}.pep.xml",
    log:
        o="log/{db}/xtandem2pepxml/{ds}/{sdb}/{xp}.out",
        e="log/{db}/xtandem2pepxml/{ds}/{sdb}/{xp}.err"
    params:
    threads: 1
    resources:
        mem = 4000
    shell:'''
        pepxmltk.py {input.txml} {output.pepxml}   
    '''

rule xtandem_group:
    '''
        Check that all xtandem output files are ready and write a flag
    '''
    input:
        ["out/{db}/xtandem/{ds}/{xp}.t.xml".format(db=db,xp=os.path.splitext(os.path.basename(fname))[0],ds=ds) for ds in config["datasets"] if config["software"]["xtandem"]["enabled"] and config["datasets"][ds]["enabled"] for db in config["active_dbs"] for fname in glob.glob("res/data/raw/{ds}/*.mzML".format(ds=ds))],
    output:
        touch("out/{db}/xtandem/{ds}.done")
