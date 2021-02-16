from tstk.io import parsefastx
from contextlib import contextmanager
import csv
import functools
import gzip
import hashlib


@contextmanager
def open_as_text(file,mode="r"):
    if mode not in {"a","r","w","x"}:
        raise ValueError("unknown file mode: '{}'".format(mode))
    text_mode = mode + "t"
    if str(file).lower().endswith(".gz"):
        opener = functools.partial(gzip.open,mode=text_mode)
    else:
        opener = functools.partial(open,mode=text_mode)
    file_obj = None
    try:
        file_obj = opener(file)
        yield file_obj
    finally:
        if file_obj is not None:
            file_obj.close()


seen = set()

counter = 0

with open(snakemake.output.db,"w") as ofh, open(snakemake.output.meta_file,"w") as mfh:
    headings = ("db_name","db_seq_id","seq_md5")
    meta_writer = csv.DictWriter(mfh,headings,dialect="excel-tab")
    meta_writer.writeheader()
    with open(snakemake.log.o,"w") as lfh:
        for db,fpath in zip(snakemake.params.dbs,snakemake.input.paths):
            with open_as_text(fpath) as fh:
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
