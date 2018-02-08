import os,glob
from tstk.io import parsefastx
import xml.etree.ElementTree as ET

configfile: "config.yaml"

rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        ["out/{db}/comet/{ds}/{xp}.pep.xml".format(db=db,xp=os.path.splitext(os.path.basename(fname))[0],ds=ds) for ds in config["datasets"] for db in config["active_dbs"] for fname in glob.glob("res/data/raw/{ds}/*.mzML".format(ds=ds))],
        ["out/{db}/xtandem/{ds}/{xp}.xml".format(db=db,xp=os.path.splitext(os.path.basename(fname))[0],ds=ds) for ds in config["datasets"] if config["datasets"][ds]["enabled"] for db in config["active_dbs"] for fname in glob.glob("res/data/raw/{ds}/*.mzML".format(ds=ds))],

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
        raw="res/data/raw/{ds}/{xp}.mzML",
        db="out/{db}/db/db_and_decoys.fasta",
    output:
        xml="out/{db}/comet/{ds}/{xp}.pep.xml"
    log:
        o="log/{db}/comet/{ds}/{xp}.out",
        e="log/{db}/comet/{ds}/{xp}.err"
    params:
        bin=config["software"]["comet"]["bin"],
        params=config["software"]["comet"]["params"],
        basename=lambda wc: "out/{db}/comet/{ds}/{xp}".format(db=wc.db,xp=wc.xp,ds=wc.ds)
    threads: config["software"]["comet"]["threads"]
    resources:
        mem = 8000
    shell:"""
        export NSLOTS={threads}
        {params.bin} -P{params.params} -D{input.db} -N{params.basename} {input.raw} > {log.o} 2> {log.e}
    """

rule xtandem_taxonomy:
    input:
        db="out/{db}/db/db_and_decoys.fasta"
    output:
        tax="out/{db}/db/xtandem.taxonomy.xml"
    run:
        root = ET.Element("biomla", label="x! taxon-to-file matching list")
        taxon = ET.SubElement(root,"taxon", label=config["dbs"][wildcards.db]["tax"])
        f = ET.SubElement(taxon,"file", format="peptide", URL=config["dbs"][wildcards.db]["path"])

        tree = ET.ElementTree(root)
        tree.write(output.tax, xml_declaration=True)

rule xtandem:
    '''
        Run xtandem on the files against the "decoyed" database.
    '''
    input:
        raw="res/data/raw/{ds}/{xp}.mgf",
        tax="out/{db}/db/xtandem.taxonomy.xml"
    output:
        xml="out/{db}/xtandem/{ds}/{xp}.xml",
        conf="out/{db}/xtandem/{ds}/{xp}.conf.xml"
    log:
        o="log/{db}/xtandem/{ds}/{xp}.out",
        e="log/{db}/xtandem/{ds}/{xp}.err"
    params:
        bin=config["software"]["xtandem"]["bin"],
        params=config["software"]["xtandem"]["params"],
        tax=config["software"]["xtandem"]["taxonomy"],
    threads: config["software"]["xtandem"]["threads"]
    resources:
        mem = 8000
    run:
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
                n.text = input.raw
            elif l == "output, path":
                n.text = output.xml

        tree.write(output.conf, xml_declaration=True)

        shell("{params.bin} {output.conf} > {log.o} 2> {log.e}")
