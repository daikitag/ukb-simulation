from utils import obtain_msprime_ratemap
import tskit
import demes
import msprime
import numpy as np
import pandas as pd
import pyslim

from snakemake.script import snakemake as snk

def subset_tree_seq(ts, selected_individuals):
    selected_nodes = np.array([], dtype=int)
    for individual in selected_individuals:
        selected_nodes = np.concatenate((selected_nodes, ts.individual(individual).nodes))

    return ts.simplify(selected_nodes, filter_individuals=False)


def convert_allele(ts):
    ts = pyslim.generate_nucleotides(ts)
    ts = pyslim.convert_alleles(ts)

    return ts

def main():
    ts = tskit.load(snk.input.ts)

    chromosome = int(snk.params.chromosome)
    arm = snk.params.arm
    recapitate_seed = int(snk.params.recapitate_seed)

    recapitate_seed *= chromosome
    recapitate_seed += 1 if arm == "p" else 0

    recombination_map, left_position = obtain_msprime_ratemap(
        recombination_map_file=snk.input.recombination_map_file,
        position_file=snk.input.position_file,
        chromosome=chromosome+arm,
    )

    demography = demes.load(snk.input.demography)

    ts = msprime.sim_ancestry(
        initial_state=ts,
        demography=demography,
        recombination_rate=recombination_map,
        random_seed=recapitate_seed,
    )
    ts = ts.simplify()

    ts = ts.shift(left_position)

    individual_id_df = pd.read_csv(snk.input.individual_id)

    ts = subset_tree_seq(ts, individual_id_df["individual_id"])

    ts = convert_allele(ts)

    ts.dump(snk.output.ts)

if __name__=="__main__":
    main()
