import collections
import sys

import numpy as np
import pandas as pd
import tskit
import tszip


def count_site_alleles(ts, tree, site):
    """Obtain collections Counter object of ancestral state and number of samples
    from the input site.
    """
    counts = collections.Counter({site.ancestral_state: ts.num_samples})
    for m in site.mutations:
        current_state = site.ancestral_state
        if m.parent != tskit.NULL:
            current_state = ts.mutation(m.parent).derived_state
        # Silent mutations do nothing
        if current_state != m.derived_state:
            num_samples = tree.num_samples(m.node)
            counts[m.derived_state] += num_samples
            counts[current_state] -= num_samples
    return counts


def obtain_maf(ts):
    maf_count = []

    tree = tskit.Tree(ts)

    for i in range(ts.num_sites):
        site = ts.site(i)
        tree.seek(site.position)
        counts = count_site_alleles(ts, tree, site)
        # counts is a Counter object from collections
        max_allele_count = counts.most_common(1)[0][1]
        freq = max_allele_count / ts.num_samples

        maf_count.append(1 - freq)

    return maf_count


def subset_tree_seq(ts, selected_individuals):
    selected_nodes = np.array([], dtype=int)
    for individual in selected_individuals:
        selected_nodes = np.concatenate(
            (selected_nodes, ts.individual(individual).nodes)
        )

    subset_ts = ts.simplify(selected_nodes)

    return subset_ts


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    ts = tszip.load(snakemake.input.ts)
    individual_id_df = pd.read_csv(snakemake.input.individual_id)

    chromosome = int(snakemake.params.chromosome)
    arm = snakemake.params.arm

    seed = int(snakemake.params.slim_seed)

    # This is used to set the seed for each chromosome and arm as a different
    # interger
    seed *= chromosome
    seed += 1 if arm == "p" else 0

    rng = np.random.default_rng(seed=seed)

    for pop in ["CEU", "CHB", "JPT", "YRI"]:
        pop_df = individual_id_df[individual_id_df.population == pop]
        pop_number = int(snakemake.params[f"{pop.lower()}_number"])
        pop_individual_id = rng.choice(
            pop_df.individual_id, size=pop_number, replace=False
        )
        pop_ts = subset_tree_seq(ts, pop_individual_id)

        pop_maf = obtain_maf(pop_ts)
        pop_maf_df = pd.DataFrame({"MAF": pop_maf})
        pop_maf_df.to_csv(snakemake.output[f"{pop.lower()}_maf"])


if __name__ == "__main__":
    main()
