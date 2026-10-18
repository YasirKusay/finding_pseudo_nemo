import warnings

print("\t".join(["CHROM", "POS", "REF", "ALT", "FILTER", "AC", "AN", "AF"]))

with open("ikbkg_ikbkgp1_gnomad_entires_VCF.tsv", "r") as f:
    for line in f:
        curr = line.strip().split("\t")
        chr = curr[0]
        pos = curr[1]
        ref = curr[3]
        alt = curr[4]
        filter = curr[6]

        # AC_joint=0;AN_joint=271642;AF_joint=0;
        INFO = curr[7].split(";")
        INFO_dict = {}

        for item in INFO:
            info_field = item.split("=")
            if info_field[0].startswith("A"):
                INFO_dict[info_field[0]] = info_field[1]

        try:
            ac = INFO_dict["AC_joint"]
            an = INFO_dict["AN_joint"]
            af = INFO_dict["AF_joint"]

            print("\t".join([chr,pos,ref,alt,filter,ac,an,af]))
        except:
            warnings.warn("Skipping: " + line)

