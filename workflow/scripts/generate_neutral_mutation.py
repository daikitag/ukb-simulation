import tskit
import msprime
import tszip

from snakemake.script import snakemake as snk

def main():
    ts = tskit.load(snk.input.ts)

    ts = msprime.sim_mutations(
        ts,
        rate=float(snk.params.neutral_mu),
        keep=True,
        random_seed=int(snk.params.neutral_seed),
    )

    tszip.compress(ts, snk.output.ts)

if __name__ == '__main__':
    main()