import sys

import msprime
import tskit
import tszip


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    ts = tskit.load(snakemake.input.ts)

    ts = msprime.sim_mutations(
        ts,
        rate=float(snakemake.params.neutral_mu),
        keep=True,
        random_seed=int(snakemake.params.neutral_seed),
    )

    tszip.compress(ts, snakemake.output.ts)


if __name__ == "__main__":
    main()
