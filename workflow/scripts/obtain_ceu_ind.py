import sys

import numpy as np
import pandas as pd


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)
    pc1 = int(snakemake.params.pc1)
    num_ceu_gwas = int(snakemake.params.num_ceu_gwas)
    pcs_df = pd.read_csv(snakemake.input.pcs, delimiter="\t")
    ceu_df = pcs_df[pcs_df["IID"].str.contains("CEU")]
    ceu_df = ceu_df[ceu_df["PC1"] < pc1]

    chromosome = int(snakemake.params.chromosome)
    arm = snakemake.params.arm

    seed = int(snakemake.params.individual_seed)

    # This is used to set the seed for each chromosome and arm as a different
    # interger
    seed *= chromosome
    seed += 1 if arm == "p" else 0

    rng = np.random.default_rng(seed=seed)

    individual_id = rng.choice(ceu_df["IID"], size=num_ceu_gwas, replace=False)
    individual_id_df = pd.DataFrame({"#FID": individual_id, "IID": individual_id})
    individual_id_df.to_csv(snakemake.output.individual_id, sep="\t", index=False)


if __name__ == "__main__":
    main()
