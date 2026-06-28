from utils import obtain_msprime_ratemap
import tskit
import demes
import msprime
import numpy as np
import pandas as pd
import pyslim

from snakemake.script import snakemake as snk

def obtain_population_id(ts):
    population_ids = {}

    for population in ts.populations():
        name = population.metadata["name"]
        if name in ["YRI", "CEU", "CHB", "JPT"]:
            population_ids[name] = population.id

    return population_ids

def obtain_individual_df(ts, yri_number, ceu_number, chb_number, jpt_number, rng):
    population_ids = obtain_population_id(ts)

    selected_by_pop = {
        "YRI": np.sort(rng.choice(np.unique(ts.nodes_individual[ts.samples(population_ids["YRI"])]), yri_number, replace=False)),
        "CEU": np.sort(rng.choice(np.unique(ts.nodes_individual[ts.samples(population_ids["CEU"])]), ceu_number, replace=False)),
        "CHB": np.sort(rng.choice(np.unique(ts.nodes_individual[ts.samples(population_ids["CHB"])]), chb_number, replace=False)),
        "JPT": np.sort(rng.choice(np.unique(ts.nodes_individual[ts.samples(population_ids["JPT"])]), jpt_number, replace=False)),
    }

    return pd.DataFrame([
        {"individual_id": ind_id, "plink_id": f"{pop}_{ind_id}", "population": pop}
        for pop, inds in selected_by_pop.items()
        for ind_id in inds
    ])

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

    # Subset tree sequence by individual IDs
    rng = np.random.default_rng(seed=int(snk.params.individual_seed))

    individual_id_df = obtain_individual_df(
        ts = ts,
        ceu_number = int(snk.params.ceu_number),
        yri_number = int(snk.params.yri_number),
        chb_number = int(snk.params.chb_number),
        jpt_number = int(snk.params.jpt_number),
        rng = rng
    )

    individual_id_df.to_csv(snk.output.individual_id, index=False)

    ts = subset_tree_seq(ts, individual_id_df["individual_id"])

    ts = convert_allele(ts)

    ts.dump(snk.output.ts)

if __name__=="__main__":
    main()
