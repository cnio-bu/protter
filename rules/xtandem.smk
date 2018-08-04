rule xtandem_taxonomy:
    input:
        db="out/{db}/db/{td}.fasta"
    output:
        tax="out/{db}/xtandem/taxonomy/{td}.xtandem.taxonomy.xml"
    conda:
        "../envs/protter.yaml"
    run:
        import xml.etree.ElementTree as ET

        root = ET.Element("biomla", label="x! taxon-to-file matching list")
        taxon = ET.SubElement(root,"taxon", label=config["dbs"][wildcards.db]["tax"])
        f = ET.SubElement(taxon,"file", format="peptide", URL=input.db)

        tree = ET.ElementTree(root)
        tree.write(output.tax, xml_declaration=True)

rule xtandem:
    '''
        Run xtandem on the files.
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
    conda:
        "../envs/protter.yaml"
    run:
        import xml.etree.ElementTree as ET

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
        pepxmltk.py {input.txml} {output.pepxml} > {log.o} 2> {log.e}  
    '''
