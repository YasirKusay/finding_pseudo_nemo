# Usage: python3 CNV-Z.py sample_depth_file
# Supply sample depth file to this command
# Easiest way is to run main_workflow.sh, then get the file inside:
# ${sample}/preprocessing/${sample}.remapped_depth.tsv

from pathlib import Path
import sys

import numpy as np
import pandas as pd

HOM_REGION_BED = "refs/beds/ikbkg_ikbkgp1_COMPLETE_HOMOLOGOUS_REGIONS.bed"

def main():
    sample_depth_file = sys.argv[1]
    sample = Path(sample_depth_file).name.split(".")[0]

    files = [sample_depth_file]
    for i in range(1, 55):
        files.append("cnv_z_depths/control_" + str(i) + ".remapped_depth.tsv")

    # ---- create a DataFrame for each file above and push it into Fx ----
    Fx = [pd.read_csv(f, sep="\t") for f in files]

    # ---- read the two homologous-region lines from the bed file ----
    with open(HOM_REGION_BED, "r") as f:
        line1 = f.readline().strip().split("\t")
        line2 = f.readline().strip().split("\t")

    # add 1 as bed is 0-based (comment carried over from original; see note below)
    HOM_1_START = int(line1[1]) + 1
    HOM_1_END = int(line1[2])

    HOM_2_START = int(line2[1]) + 1
    HOM_2_END = int(line2[2])

    # ---- excluding the 35kb homologous regions ----
    # chose this region to exclude, as we want a stable CN=2 region to normalise against
    # allowing cleaner comparisons between other samples
    # as long as normalisation method is consistent, it will be fine
    def normgene(chr_, pos):
        ichr = chr_ == "chrX"
        ipos = not ((HOM_1_START <= pos <= HOM_1_END) or (HOM_2_START <= pos <= HOM_2_END))
        return ichr and ipos

    # Need to adjust this such that we are only normalising against our
    # region of interest
    for i in range(len(Fx)):
        mask = [normgene(c, p) for c, p in zip(Fx[i]["Chr"], Fx[i]["Position"])]
        denom = Fx[i].loc[mask, "Depth"].sum()
        Fx[i]["prop"] = Fx[i]["Depth"] / denom

    # ---- build the M matrix (rows = positions, cols = samples) ----
    M = np.empty((len(Fx[0]), len(Fx)), dtype=np.float64)
    for i in range(len(Fx)):
        M[:, i] = Fx[i]["prop"].to_numpy()

    # ---- row-wise mean across samples ----, axis=1 means rows
    row_mean = M.mean(axis=1)
    for i in range(len(Fx)):
        Fx[i]["mean"] = row_mean

    # ---- row-wise std across samples (ddof=1 to match Julia's Statistics.std) ----
    row_std = M.std(axis=1, ddof=1)
    for i in range(len(Fx)):
        Fx[i]["std"] = row_std

    # ---- expected depth ----
    for i in range(len(Fx)):
        Fx[i]["expDepth"] = (Fx[i]["Depth"] / Fx[i]["prop"] * Fx[i]["mean"])

    # ---- z-score ----
    for i in range(len(Fx)):
        Fx[i]["zscore"] = ((Fx[i]["prop"] - Fx[i]["mean"]) / Fx[i]["std"])

    # ---- copy number ----
    for i in range(len(Fx)):
        Fx[i]["copynumber"] = 2 * Fx[i]["prop"] / Fx[i]["mean"]

    # ---- write output file ----
    # index 0 is our sample of interest
    Fx[0].to_csv(f"Fx-{sample}.csv", index=False)

if __name__ == "__main__":
    main()