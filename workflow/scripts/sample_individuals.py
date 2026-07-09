import sys

import numpy as np
import pandas as pd
import tskit


def obtain_population_id(ts):
    """Obtain population IDs for four populations.

    This function obtains population IDs for CEU, YRI, CHB and JPT from the
    input tree sequence data.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Input tree sequence data.

    Returns
    -------
    population_ids : dict
        Dictionary with population as a key and population ID as a value.
    """
    population_ids = {}

    for population in ts.populations():
        name = population.metadata["name"]
        if name in ["YRI", "CEU", "CHB", "JPT"]:
            population_ids[name] = population.id

    return population_ids


def obtain_individual_df(ts, yri_number, ceu_number, chb_number, jpt_number, rng):
    """Obtain individual ID dataframe.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Input tree sequence data.
    yri_number : int
        Number of YRI individuals to be selected.
    ceu_number : int
        Number of CEU individuals to be selected.
    chb_number : int
        Number of CHB individuals to be selected.
    jpt_number : int
        Number of JPT individuals to be selected.
    rng : numpy.random.Generator
        Random generator that will be used to select individuals.

    Returns
    -------
    individual_df : pandas.DataFrame
        Dataframe with individual IDs that are selected from the tree sequence data.
        There are three columns in this dataframe:
        - `individual_id` : This is the individual ID in the tree sequence data.
        - `plink_id` : This is named as {population}_{individual ID}. For example,
          if an individual 1 from CEU is being selected, the PLINK ID will be
          CEU_1. This ID will be used to generate PLINK files, and it is to make
          sure that we can easily classify individual's population from their name.
        - `population` : Population of an individual.
    """
    population_ids = obtain_population_id(ts)

    selected_by_pop = {
        "YRI": np.sort(
            rng.choice(
                np.unique(ts.nodes_individual[ts.samples(population_ids["YRI"])]),
                yri_number,
                replace=False,
            )
        ),
        "CEU": np.sort(
            rng.choice(
                np.unique(ts.nodes_individual[ts.samples(population_ids["CEU"])]),
                ceu_number,
                replace=False,
            )
        ),
        "CHB": np.sort(
            rng.choice(
                np.unique(ts.nodes_individual[ts.samples(population_ids["CHB"])]),
                chb_number,
                replace=False,
            )
        ),
        "JPT": np.sort(
            rng.choice(
                np.unique(ts.nodes_individual[ts.samples(population_ids["JPT"])]),
                jpt_number,
                replace=False,
            )
        ),
    }

    individual_df = pd.DataFrame(
        [
            {"individual_id": ind_id, "plink_id": f"{pop}_{ind_id}", "population": pop}
            for pop, inds in selected_by_pop.items()
            for ind_id in inds
        ]
    )

    return individual_df


def main():
    """Obtain individual ID dataframe.

    This function will extract a pandas Dataframe with randomly selected
    individuals from CEU, YRI, CHB, and JPT populations.
    """
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    ts = tskit.load(snakemake.input.ts)

    rng = np.random.default_rng(seed=int(snakemake.params.individual_seed))

    individual_id_df = obtain_individual_df(
        ts=ts,
        ceu_number=int(snakemake.params.ceu_number),
        yri_number=int(snakemake.params.yri_number),
        chb_number=int(snakemake.params.chb_number),
        jpt_number=int(snakemake.params.jpt_number),
        rng=rng,
    )

    individual_id_df.to_csv(snakemake.output.individual_id, index=False)


if __name__ == "__main__":
    main()
