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
    conda:
        "../envs/protter.yaml"
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
    conda:
        "../envs/protter.yaml"
    run:
        import numpy as np
        import pandas as pd
        from statsmodels.stats.multitest import fdrcorrection_twostage

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

