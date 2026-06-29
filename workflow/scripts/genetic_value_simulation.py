import numpy as np
import pandas as pd
import tskit
import tstrait

import sys

def obtain_mutation_df(ts, selection_scaling):
    """Obtain mutation dataframe.

    Obtain mutation dataframe from a tree sequence that is simulated in slim
    We will multiply the final selection coefficient by -2*selection_scaling, due to
    how we are modeling selection coefficient in SLiM simulation.

    Parameters
    ----------
    ts : tskit.TreeSequence
        The tree sequence data that will be used to compute mutation dataframe.
    selection_scaling : float
        Selection scaling parameter in SLiM simulation.

    Returns
    -------
    mutation_df : pandas.DataFrame
        Mutation dataframe from the input tree sequence data.

    Notes
    -----
    Only causal mutations are required in the input tree sequence data, as we extract
    information regarding selection coefficient and causal allele. Neutral mutations
    are not required. Be sure to input tree sequence data that has the alleles in the
    desired format (i.e., nucleotides and not SLiM allele if the final output uses
    nucleotides), as the alleles inside this dataframe will be matched to the final output.

    The output mutation dataframe includes the following columns. These information are
    extracted for all mutations in the input tree sequence data:

    - `site_id`: Site ID from the tree sequence data.
    - `site_position`: Location of the mutation in the input tree sequence data.
    - `selection_coeff`: Selection coefficient that is simulated in SLiM underdominant
      selection simulation.
    - `causal_allele`: Derived state of the mutation. This will be used as a causal
      alllele in phenotypic simulation.
    - `trait_id`: This will be set as 0, because we are only simulating a single trait.
      If you are interested in simulating multiple traits, set the seed of genetic
      value simulation as a different number.
    """
    mutation_list = []
    for i in range(ts.num_mutations):
        mutation = ts.mutation(i)

        mutation_list.append(
            {
                "site_id": ts.site(mutation.site).id,
                "site_position": ts.site(mutation.site).position,
                "selection_coeff": mutation.metadata["mutation_list"][0]["selection_coeff"],
                "causal_allele": mutation.derived_state,
                "trait_id": 0,
            }
        )

    mutation_df = pd.DataFrame(mutation_list)
    mutation_df["selection_coeff"] = mutation_df["selection_coeff"]*selection_scaling*(-2)

    return mutation_df

def sim_tstrait_pleiotropy(ts, selection_scaling, proportion, seed):
    """Obtain mutation dataframe with simulated effect sizes and take a subset of them.

    Parameters
    ----------
    ts : tskit.TreeSequence
        The tree sequence data that will be used to simulate effect sizes and obtain
        mutation information.
    selection_scaling : float
        Selection scaling parameter in SLiM simulation.
    proportion : float
        Proportion of the mutation dataframe that will be included in the final
        output. This should not be larger than `1`.
    seed : int
        Seed that will be used to simulate effect sizes.

    Returns
    -------
    subset_mutation_df : pandas.DataFrame
        A subset of the mutation dataframe from the input tree sequence data.
    genetic_df : pandas.DataFrame
        A dataframe that includes information regarding genetic values of all
        individuals in the input tree sequence data.

    Notes
    -----
    This function is comprised of four steps, (1) obtaining information regarding
    mutations, (2) simulating effect sizes by using the selection coefficient of
    each mutation, (3) taking a subset of that mutation dataframe, and (4)
    compuate genetic value of individuals by using the mutation dataframe.
    
    In step (1), :func:`obtain_mutation_df` function is used to obtain a mutation
    dataframe by extracting information regarding how mutation dataframe is obtained
    from the input tree sequence data.

    In step (2), effect sizes are simulated from a normal distribution with mean `0`
    and variance `selection_coefficient`. The influence of heritability and target size
    of the trait are being modeled elsewhere, so we simply let the variance be
    `selection_coefficient`.

    In step (3), a proportion of the mutation dataframe is selected based on the
    `proportion` parameter. While all mutation in the tree sequence data were used in
    SLiM's underdominant selection simulation, a subset of them will be selected to
    conduct quantitative trait simulation to account for the mutation target size of
    the trait. Random sampling without replacement is used to extract a proportion of
    the resulting mutation dataframe.

    In step (4), :func:`tstrait.genetic_value` is used to compute the genetic values of all
    individuals in the input tree sequence data.

    The input seed is directly used to simulate effect sizes, but a seed of
    `seed + 1` is used to take a subset of the mutation dataframe.

    The output mutation dataframe includes the following columns. Once these information
    are extracted for all mutations, a subset of them are selected in the final output.
    This mutation dataframe is ordered by `site_id`, as it is a requirement in `tstrait`
    genetic value computation. This mutation dataframe is used to compute the genetic value
    of individuals.

    - `site_id`: Site ID from the tree sequence data.
    - `site_position`: Location of the mutation in the input tree sequence data.
    - `selection_coeff`: Selection coefficient that is simulated in SLiM underdominant
      selection simulation.
    - `causal_allele`: Derived state of the mutation. This will be used as a causal
      alllele in phenotypic simulation.
    - `trait_id`: This will be set as 0, because we are only simulating a single trait.
      If you are interested in simulating multiple traits, set the seed of genetic
      value simulation as a different number.
    - `effect_size`: Simulated effect size.

    The genetic value dataframe is simply obtained by using
    :func:`tstrait.genetic_value` function. Please refer to the :func:`tstrait.genetic_value`
    documentation to find out the details of the output dataframe.
    """
    rng = np.random.default_rng(seed=seed)
    mutation_df = obtain_mutation_df(ts, selection_scaling)
    mutation_df["effect_size"] = mutation_df.apply(lambda row: rng.normal(loc = 0, scale = np.sqrt(row["selection_coeff"])), axis=1)

    if proportion > 1:
        raise ValueError("Proportion must not be greater than 1.")

    subset_mutation_df = mutation_df.sample(
        frac=proportion, replace=False, random_state=seed+1
    )
    subset_mutation_df = subset_mutation_df.sort_values(by="site_id")

    genetic_df = tstrait.genetic_value(ts, subset_mutation_df)

    return subset_mutation_df, genetic_df


def main():
    """Run the genetic value simulation workflow.

    This function will extract two `pandas` DataFrame, mutation dataframe
    and genetic value dataframe, from the input tree sequence data. The
    mutation dataframe includes information regarding mutations that are used in
    quantitative trait simulation, and the genetic value dataframe includes
    information regarding all individuals in the tree sequence data.
    """
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    ts = tskit.load(snakemake.input.ts)

    subset_mutation_df, genetic_df = sim_tstrait_pleiotropy(
        ts = ts,
        selection_scaling = float(snakemake.params.selection_scaling),
        seed = int(snakemake.params.genetic_seed),
        proportion = float(snakemake.params.proportion),
    )

    subset_mutation_df.to_csv(snakemake.output.mutation_df, index=False)

    genetic_df.to_csv(snakemake.output.genetic_df, index=False)

if __name__ == '__main__':
    main()