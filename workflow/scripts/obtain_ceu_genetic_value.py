import sys

import pandas as pd


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)
    genetic_df = pd.read_csv(snakemake.input.genetic_df)
    individual_df = pd.read_csv(snakemake.input.individual_id, sep="\t")

    ts_ind_id = [int(x.split("_")[1]) for x in individual_df["IID"]]
    genetic_df = genetic_df[genetic_df["individual_id"].isin(ts_ind_id)]

    genetic_df.to_csv(snakemake.output.genetic_df, index=False)


if __name__ == "__main__":
    main()
