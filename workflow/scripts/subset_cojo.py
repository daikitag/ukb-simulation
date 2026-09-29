import sys

import pandas as pd


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)
    cojo_df = pd.read_csv(snakemake.input.cojo_result, sep="\t")
    cojo_df = cojo_df[cojo_df.pJ < 5e-8]
    cojo_df = cojo_df[["SNP", "A1", "bJ"]]

    cojo_df.to_csv(snakemake.output.cojo_df, sep="\t", index=False)


if __name__ == "__main__":
    main()
