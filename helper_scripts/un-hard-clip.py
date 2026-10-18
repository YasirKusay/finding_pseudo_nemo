import sys, gzip, pysam

COMP = str.maketrans("ACGTNacgtnRYKMBVDHrykmbvdh", "TGCANtgcanYRMKVBHDyrmkvbhd")
revcomp = lambda s: s.translate(COMP)[::-1]

fq = {}
for i, path in enumerate(sys.argv[2:5]):
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as fh:
        for j, line in enumerate(fh):
            if j % 4 == 0:
                name = line[1:].strip()
                if name not in fq.keys():
                    fq[name] = []
            elif j % 4 == 1:
                fq[name].append(line.strip())
 
fin = pysam.AlignmentFile(sys.argv[1])
with pysam.AlignmentFile("-", "w", template=fin) as fout:
    for r in fin:
        if r.is_supplementary and "H" in (r.cigarstring or ""):

            if r.is_read2:
                seq = fq[r.query_name][1]
            else:
                seq = fq[r.query_name][0]
            if r.is_reverse:
                seq = revcomp(seq)
            assert r.query_sequence in seq

            r.cigarstring = r.cigarstring.replace("H", "S")
            r.query_sequence = seq

        fout.write(r)