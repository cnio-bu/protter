import os,glob,math
import numpy as np
from tstk.io import parsefastx
import xml.etree.ElementTree as ET
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection_twostage

configfile: "config.yaml"

def get_samples(ds,grouping):
    fpath = "res/data/prot/{ds}".format(ds=ds)
    #get the grouping statement from the config file and create a "gf" function with it
    exec('gf = lambda x: {}'.format(config["datasets"][ds]["groupings"][grouping]), globals())
    groups = {}
    for fname in glob.glob("{}/*.mzML".format(fpath)):
        f = os.path.splitext(os.path.basename(fname))[0]
        g = gf(f)
        try:
            groups[g].append(f)
        except KeyError:
            groups[g] = [f]
    return groups

def input_crux_percolator():
    for software in config["software"]["search"]:
        if config["software"]["search"][software]["enabled"]:
            for ds in config["datasets"]:
                if config["datasets"][ds]["enabled"]:
                    for grouping in config["datasets"][ds]["groupings"]:
                        for group in get_samples(ds,grouping):
                            for db in config["dbs"]:
                                if config["dbs"][db]["enabled"] and db in config["datasets"][ds]["dbs"]:
                                    for em in ["tdc","mix-max"]:
                                        yield "out/{db}/percolator/{sw}/{ds}/{em}/{grouping}/{group}/percolator.target.peptides.txt".format(grouping=grouping,group=group,sw=software,db=db,ds=ds,em=em)

def input_add_confidence():
    for software in config["software"]["search"]:
        if config["software"]["search"][software]["enabled"]:
            for ds in config["datasets"]:
                if config["datasets"][ds]["enabled"]:
                    for db in config["dbs"]:
                        if "add_confidence" in config["dbs"][db]:
                            if config["dbs"][db]["enabled"] and db in config["datasets"][ds]["dbs"]:
                                    for sample in get_samples(ds,'single'):
                                        yield "out/{db}/add_confidence/{sw}/{ds}/{sample}.tsv".format(sw=software,db=db,ds=ds,sample=sample)

rule all:
    '''
        Main rule, which requires as input the final output of the workflow.
    '''
    input:
        input_crux_percolator(),
        #input_add_confidence()

rule process_db:
    input:
        dbs=lambda wc: config["dbs"][wc.db]["paths"]
    output:
        db="out/{db}/db/target.fasta"
    run:
        with open(output.db,"w") as ofh:
            for fpath in input.dbs:
                with open(fpath) as fh:
                    for r in parsefastx(fh):
                        ofh.write(">{}\n{}\n".format(r[0],r[1].replace("L","I")))

rule add_decoys:
    '''
        Use decoPYrat to generate decoy sequences.
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
        python bin/decoyPYrat.py {input.db} -t {params.tmp}  -o {output.dec} --decoy_prefix decoy -k > {log.o} 2> {log.e}
        cat {input.db} {output.dec} >> {output.dbd} 2>> {log.e}
    """

rule comet:
    '''
        Run comet on the files.
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

rule psm_convert:
    input:
        db="out/{db}/db/{td}.fasta",
        psm="out/{db}/{sw}/{ds}/{sample}.{td}.pep.xml"
    output:
        d=temp("out/{db}/{sw}/{ds}/{sample}.{td}.psm-convert"),
        f="out/{db}/{sw}/{ds}/{sample}.{td}.tsv"
    log:
        o="log/{db}/psm-convert/{ds}/{sample}.{td}.out",
        e="log/{db}/psm-convert/{ds}/{sample}.{td}.err"
    params:
    threads: 1
    resources:
        mem = 1000
    shell:"""
        bin/crux/crux psm-convert --overwrite T --protein-database {input.db} --output-dir {output.d} {input.psm} tsv > {log.o} 2> {log.e}
        mv {output.d}/psm-convert.txt {output.f}
        rm {output.d}/*
    """
    
def dfilter(df,field,strings,neg=True):
    '''Function that filters a dataframe based on a field and a number of strings'''
    if isinstance(strings, str):
        strings = [strings]
    for string in strings:
        if neg:
            df = df[~df[field].str.contains(string)]
        else:
            df = df[df[field].str.contains(string)]
    return df

rule add_confidence:
    input:
        t="out/{db}/{sw}/{ds}/{sample}.target.tsv",
        d="out/{db}/{sw}/{ds}/{sample}.decoy.tsv"
    output:
        f="out/{db}/add_confidence/{sw}/{ds}/{sample}.tsv"
    params:
    threads: 1
    benchmark: "benchmark/{db}/add_confidence/{sw}/{ds}/{sample}.txt"
    resources:
        mem = 32000
    run:
        filters  = config["dbs"][wildcards.db]["add_confidence"]

        #run filter separately to save memory
        dft = pd.read_table(input.t)
        for f in filters:
            neg = f["neg"]
            col = f["column"]
            strings = f["strings"]
            dft = dfilter(dft,col,strings,neg)

        dfd = pd.read_table(input.d)
        for f in filters:
            neg = f["neg"]
            col = f["column"]
            strings = f["strings"]
            dfd = dfilter(dfd,col,strings,neg)

        df = dft.append(dfd,ignore_index = True)

        #corrected xcorr as provided by cnic
        df["R"] = df.apply(lambda x: 1.0 if x["charge"] < 3 else 1.2, axis=1)
        df["xcorr2"] = df.apply(lambda x: np.nan if x["xcorr score"] == 0 else math.log10(x["xcorr score"] / x["R"]) / math.log10(2 * x["peptide mass"] / 110), axis=1) 
        df = df.drop('R', 1)

        #FDR calculations
        df["is_decoy"] = df.apply(lambda x: 1 if "decoy" in x["protein id"] else 0, axis=1)
        df["is_target"] = (df["is_decoy"] - 1) * -1
        total_decoys = df["is_decoy"].sum()

        #FDR calculations by PSM
        for f in ["sp score","xcorr2"]:
            #define field names
            ndecoysf = "ndecoys_{}".format(f) #cumulative sum of decoys
            ntargetsf = "ntargets_{}".format(f) #cumulative sum of targets
            pvalf = "pval_{}".format(f) #p value
            fdrf = "fdr_{}".format(f) #FDR
            padjf = "padj_{}".format(f) #FDR adjusted p-value

            df = df.sort_values(f, ascending=False)

            df[ndecoysf] = df["is_decoy"].cumsum()
            df[ntargetsf] = df["is_target"].cumsum()
            df[pvalf] = df[ndecoysf] / total_decoys
            df[fdrf] = df[ndecoysf] / df[ntargetsf]

            #adjusted p_value
            df[padjf] = fdrcorrection_twostage(df[pvalf])[1]
            
        #FDR calculations by peptide
        for f in ["sp score","xcorr2"]:
            seen_pep = set() #store peptides that we've already seen
            ndecoysf = "ndecoys_{}_pep".format(f)
            ntargetsf = "ntargets_{}_pep".format(f)
            pvalf = "pval_{}_pep".format(f)
            fdrf = "fdr_{}_pep".format(f)
            padjf = "padj_{}_pep".format(f)

            df = df.sort_values(f, ascending=False)
            
            n = {"tgt":[],"dec":[]} # lists to store cumulative sum of "novel" peptides
            
            for index, row in df.iterrows():
                if row["is_target"] == 1:
                    dest = "tgt"
                elif row["is_decoy"] == 1:
                    dest = "dec"
                else:
                    raise ValueError("Row is neither target nor decoy {}".format(row))
                
                #initialise cumsums stores in first row
                #TODO: is there a better way to do this?
                if not n[dest]:
                    if row["is_target"] == 1:
                        n["tgt"].append(1)
                        n["dec"].append(0)
                    elif row["is_decoy"] == 1:
                        n["tgt"].append(0)
                        n["dec"].append(1)
                    seen_pep.add(row["sequence"])
                    continue
                    
                #keep the same cumsum as before for the current row
                n["tgt"].append(n["tgt"][-1])
                n["dec"].append(n["dec"][-1])

                #increase the corresponding cumsum by one if novel peptide
                if row["sequence"] not in seen_pep:
                    seen_pep.add(row["sequence"])
                    n[dest].append(n[dest].pop()+1)
                        
            df[ndecoysf] = pd.Series(n["dec"]).values
            df[ntargetsf] = pd.Series(n["tgt"]).values
            df[pvalf] = df[ndecoysf] / total_decoys
            df[fdrf] = df[ndecoysf] / df[ntargetsf]

            #adjusted p_value
            df[padjf] = fdrcorrection_twostage(df[pvalf])[1]

        df.to_csv(output.f,sep="\t")

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
        pepxmltk.py {input.txml} {output.pepxml} > {log.o} 2> {log.e}  
    '''

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
    resources:
        mem = 64000
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

rule crux_assign_confidence:
    '''
        Run the crux assign-confidence agorithm against PSMs.
    '''
    input:
        targets = lambda wc: ["out/{db}/{sw}/{ds}/{sample}.target.pep.xml".format(db=wc.db,sw=wc.sw,ds=wc.ds,sample=sample) for sample in get_samples(wc.ds,wc.grouping)[wc.group]],
        decoys = lambda wc: ["out/{db}/{sw}/{ds}/{sample}.decoy.pep.xml".format(db=wc.db,sw=wc.sw,ds=wc.ds,sample=sample) for sample in get_samples(wc.ds,wc.grouping)[wc.group]],
    output:
        d="out/{db}/crux_assign_confidence/{sw}/{ds}/{em}/{grouping}/{group}",
        f="out/{db}/crux_assign_confidence/{sw}/{ds}/{em}/{grouping}/{group}/assign-confidence.target.txt"
    log:
        o="log/{db}/crux_assign_confidence/{sw}/{ds}/{em}/{grouping}/{group}.out",
        e="log/{db}/crux_assign_confidence/{sw}/{ds}/{em}/{grouping}/{group}.err"
    threads: 3
    resources:
        mem = 8000
    shell:"""
        bin/crux/crux assign-confidence --overwrite T --estimation-method {wildcards.em} --output-dir {output.d} {input.targets} > {log.o} 2> {log.e}
    """
