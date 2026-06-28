import numpy as np
import pandas as pd
import tskit
import tstrait
import sys

from snakemake.script import snakemake as snk

def obtain_mutation_df(ts, selection_scaling):
    """
    Obtain mutation dataframe from a tree sequence that is simulated in slim
    We will multiply the final selection coefficient by -2*selection_scaling, due to
    how we are modeling selection coefficient in SLiM simulation

    If we let `s` as the selection coefficient in the stabilizing selection model and
    model fitness of individuals in the underdominance selection model as `1`, `1 + t`,
    and `1`, we can calculate that `s = -2 * t` (see computations in the paper).

    In tstrait simulation, we must need the following columns in trait dataframe to
    compute genetic values:
    - site_id : Site ID from tree sequence data
    - effect_size : This will be simulated based on selection coefficient of mutations
    - causal_allele : This can be found from derived_state in mutation
    - trait_id : This will be set to 0, because we are only simulating a single trait
    """
    mutation_list = []
    for i in range(ts.num_mutations):
        mutation = ts.mutation(i)

        mutation_list.append(
            {
                "site_id": mutation.site,
                "selection_coeff": mutation.metadata["mutation_list"][0]["selection_coeff"],
                "causal_allele": mutation.derived_state,
                "trait_id": 0,
                "site_position": ts.site(mutation.site).position,
                "site_id": ts.site(mutation.site).id,
            }
        )

    mutation_df = pd.DataFrame(mutation_list)
    mutation_df["selection_coeff"] = mutation_df["selection_coeff"]*selection_scaling*(-2)

    return mutation_df

def sim_tstrait_pleiotropy(ts, selection_scaling, seed):
    """
    ts is the tree sequence of interest
    n is the degree of pleiotropic effects
    seed is the seed that will be used in the tstrait simulation
    w2 is the stabilizing selection parameter, which will be 1 by default
    selectin_scaling is the scaling factor for the underdominance simulation model
    """
    rng = np.random.default_rng(seed=seed)
    mutation_df = obtain_mutation_df(ts, selection_scaling)
    mutation_df["effect_size"] = mutation_df.apply(lambda row: rng.normal(loc = 0, scale = np.sqrt(row["selection_coeff"])), axis=1)
    return mutation_df

def simulate_genetic_value(ts, mutation_df):
    genetic_df = tstrait.genetic_value(ts, mutation_df)

    return genetic_df

def main():
    ts = tskit.load(snk.input[0])

    mutation_df = sim_tstrait_pleiotropy(
        ts = ts,
        selection_scaling = float(snk.params.selection_scaling),
        proportion = float(snk.params.proportion),
        seed = int(snk.params.genetic_seed),
    )

    subset_mutation_df = mutation_df.sample(
        frac=float(snk.params.proportion), replace=False, random_state=int(snk.params.genetic_seed)+1
    )
    subset_mutation_df = subset_mutation_df.sort_values(by="site_id")

    subset_mutation_df.to_csv(snk.output.mutation_df, index=False)

    genetic_df = tstrait.genetic_value(ts, subset_mutation_df)

    genetic_df.to_csv(snk.output.genetic_df, index=False)

if __name__ == '__main__':
    main()