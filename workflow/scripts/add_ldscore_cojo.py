import pandas as pd


def main():
    ld_score_df = pd.read_csv(snakemake.input.ld_score, sep=" ")
    cojo_df = pd.read_csv(snakemake.input.cojo_result, sep="\t")

    cojo_df = cojo_df.merge(ld_score_df[["SNP", "ldscore"]], on="SNP", how="left")

    cojo_df.to_csv(snakemake.output.cojo_ld_score)


if __name__ == "__main__":
    main()
