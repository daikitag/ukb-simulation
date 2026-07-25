import sys

import pandas as pd


def obtain_plink_phenotype(phenotype_df, individual_id_df, phenotype_name):
    plink_phenotype_df = pd.DataFrame(
        {
            "FID": individual_id_df["#FID"].tolist(),
            "IID": individual_id_df["IIF"].tolist(),
            phenotype_name: phenotype_df["phenotype"].to_list(),
        }
    )

    return plink_phenotype_df


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)
    phenotype_df = pd.read_csv(snakemake.input.phenotype_df)
    individual_df = pd.read_csv(snakemake.input.individual_id, sep="\t")

    plink_phenotype = obtain_plink_phenotype(
        phenotype_df=phenotype_df,
        individual_id_df=individual_df,
        phenotype_name=snakemake.params.gwas_phenotype,
    )

    plink_phenotype.to_csv(snakemake.output.plink_phenotype, sep="\t", index=False)


if __name__ == "__main__":
    main()
