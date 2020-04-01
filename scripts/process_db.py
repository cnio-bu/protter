from tstk.io import parsefastx
import csv
import hashlib

seen = set()

counter = 0

with open(snakemake.output.db,"w") as ofh, open(snakemake.output.meta_file,"w") as mfh:
    headings = ("db_name","db_seq_id","seq_md5")
    meta_writer = csv.DictWriter(mfh,headings,dialect="excel-tab")
    meta_writer.writeheader()
    with open(snakemake.log.o,"w") as lfh:
        for db,fpath in zip(snakemake.params.dbs,snakemake.input.paths):
            with open(fpath) as fh:
                for r in parsefastx(fh):
                    counter += 1

                    sid = r[0]
                    seq = r[1].replace("L","I")
                    m = hashlib.md5()
                    m.update(seq.encode('utf-8'))
                    seq_md5 = m.hexdigest()
                    meta_writer.writerow({"db_name": db, "db_seq_id": sid, "seq_md5": seq_md5})
                    if seq_md5 in seen:
                        lfh.write("{} sequence ({}) is duplicated. Ignoring.\n".format(sid,db))
                    else:
                        seen.add(seq_md5)
                        sid = "seq{}#{}#".format(counter,db)
                        ofh.write(">{}\n{}\n".format(sid,seq))
