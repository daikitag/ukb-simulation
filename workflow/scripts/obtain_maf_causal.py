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
    causal_list = []

    tree = tskit.Tree(ts)

    for i in range(ts.num_sites):
        site = ts.site(i)
        tree.seek(site.position)
        counts = count_site_alleles(ts, tree, site)
        # counts is a Counter object from collections
        max_allele_count = counts.most_common(1)[0][1]
        freq = max_allele_count / ts.num_samples

        maf_count.append(1 - freq)

        if ts.site(i).mutations[0].metadata["mutation_list"] != []:
            causal_list.append(1)
        else:
            causal_list.append(0)

    return maf_count, causal_list


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

    for pop in ["CEU", "CHB", "JPT", "YRI"]:
        pop_df = individual_id_df[individual_id_df.population == pop]
        pop_ts = subset_tree_seq(ts, pop_df.individual_id)

        pop_maf, pop_causal = obtain_maf(pop_ts)
        pop_maf_df = pd.DataFrame({"MAF": pop_maf, "causal": pop_causal})
        pop_maf_df = pop_maf_df[pop_maf_df.MAF >= 0.01]
        pop_maf_df.to_csv(snakemake.output[f"{pop.lower()}_maf_causal"], index=False)


if __name__ == "__main__":
    main()
