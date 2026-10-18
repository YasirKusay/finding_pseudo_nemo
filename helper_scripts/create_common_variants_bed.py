with open("gnomad_high_freq.tsv", "r") as f:
    for line in f:
        curr = line.strip().split("_")
        chrom = curr[0]
        start_0_base = int(curr[1]) - 1
        end = int(curr[1]) + len(curr[2]) - 1
        info = "->".join(curr[2:])
        print("\t".join([chrom, str(start_0_base), str(end), info]))
