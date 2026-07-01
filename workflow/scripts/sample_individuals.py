import tskit
import numpy as np
import pandas as pd

import sys


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

    return pd.DataFrame(
        [
            {"individual_id": ind_id, "plink_id": f"{pop}_{ind_id}", "population": pop}
            for pop, inds in selected_by_pop.items()
            for ind_id in inds
        ]
    )


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    ts = tskit.load(snakemake.input.ts)

    # Subset tree sequence by individual IDs
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
