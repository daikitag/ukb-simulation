import sys

import numpy as np
import pandas as pd


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    individual_id_df = pd.read_csv(snakemake.input.individual_id)

    score_df = pd.read_csv(snakemake.input.score_df, sep="\t")
    score_df["plink_id"] = score_df["IID"]
    score_df = score_df[
        ["plink_id", "ALLELE_CT", "NAMED_ALLELE_DOSAGE_SUM", "SCORE1_AVG"]
    ]

    chromosome = int(snakemake.params.chromosome)
    arm = snakemake.params.arm

    seed = int(snakemake.params.individual_seed) * 1000

    # This is used to set the seed for each chromosome and arm as a different
    # interger
    seed *= chromosome
    seed += 1 if arm == "p" else 0

    rng = np.random.default_rng(seed=seed)

    genetic_df = pd.read_csv(snakemake.input.genetic_df)

    # Dataframe for CHB, JPT and YRI
    for pop in ["CHB", "JPT", "YRI"]:
        pop_df = individual_id_df[individual_id_df.population == pop]
        pop_number = int(snakemake.params[f"{pop.lower()}_number"])

        pop_df = pop_df.sample(n=pop_number, replace=False, random_state=rng)

        pop_df = pd.merge(
            pop_df, genetic_df[["individual_id", "genetic_value"]], on="individual_id"
        )

        pop_df = pd.merge(pop_df, score_df, on="plink_id")

        pop_df.to_csv(snakemake.output[f"{pop.lower()}_pgs"], index=False)

    # Dataframe for CEU
    # This is to remove CEU individuals that were in the original GWAS sample
    pc1 = int(snakemake.params.pc1)

    pcs_df = pd.read_csv(snakemake.input.pcs, delimiter="\t")
    ceu_pcs_df = pcs_df[pcs_df["IID"].str.contains("CEU")]
    ceu_id = ceu_pcs_df[ceu_pcs_df["PC1"] < pc1]["IID"].to_numpy()

    ceu_gwas_id = pd.read_csv(snakemake.input.ceu_gwas_ind_id, sep="\t")[
        "IID"
    ].to_numpy()

    ceu_candidates = ceu_id[~np.isin(ceu_id, ceu_gwas_id)]

    selected_ceu_iid = rng.choice(
        ceu_candidates, size=int(snakemake.params.ceu_number), replace=False
    )

    ceu_df = individual_id_df[individual_id_df["plink_id"].isin(selected_ceu_iid)]

    ceu_df = pd.merge(
        ceu_df, genetic_df[["individual_id", "genetic_value"]], on="individual_id"
    )

    ceu_df = pd.merge(ceu_df, score_df, on="plink_id")

    ceu_df.to_csv(snakemake.output["ceu_pgs"], index=False)


if __name__ == "__main__":
    main()
