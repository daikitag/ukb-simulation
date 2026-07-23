import collections
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import tskit
import tszip


def drop_mutations(tables, indexes_of_mutations_to_keep):
    """Drop mutations from mutation table.

    This function directly modifies the table collection from a tree sequence data and
    only keeps the IDs of mutations that are specified in
    `indexes_of_mutations_to_keep`. All other mutations are removed from the table
    collection.

    Parameters
    ----------
    tables : tskit.TableCollection
        Table collection from a tree sequence data.
    indexes_of_mutations_to_keep : list
        List of mutation IDs to keep in the tree sequence data.
    """
    m = len(tables.mutations)
    tables.mutations.parent = np.zeros(m, dtype=np.int32) - 1  # null the parent column
    select = np.zeros(m, dtype=bool)
    select[indexes_of_mutations_to_keep] = True
    tables.mutations.keep_rows(select)
    tables.compute_mutation_parents()


def common_mutation_id(site, state):
    """Obtain mutation IDs with a certain state.

    The input of this function is a site from a tree sequence data and a state from
    that site. The output of this function is a list of all mutation IDs with the
    derived state equal to `state`.
    """
    mutation_list = []
    for m in site.mutations:
        if m.derived_state in state:
            mutation_list.append(m.id)
    return mutation_list


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


def maf_threshold(ts, maf):
    """Subset tree sequence based on MAF.

    This function removes all sites where the MAF is less than `maf`.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Input tree sequence data.
    maf : float
        MAF threshold.

    Returns
    -------
    down_sample_ts : tskit.TreeSequence
        Tree sequence data where the MAF of all mutations in the data will be above
        `maf`.

    Notes
    -----
    For multiallelic sites, this function only retains the mutation where the derived
    state is the major allele or the second most common allele. All other mutations
    from the site will be removed, so all sites will only have two alleles, the derived
    state from the retained mutation and the ancestral state. This is to ensure that
    the output tree sequence data is biallelic.
    """
    remove_site = []
    keep_mutation = []

    tree = tskit.Tree(ts)

    for i in range(ts.num_sites):
        site = ts.site(i)
        tree.seek(site.position)
        counts = count_site_alleles(ts, tree, site)
        # counts is a Counter object from collections
        max_allele_count = counts.most_common(1)[0][1]
        freq = max_allele_count / ts.num_samples
        if freq > (1 - maf):
            remove_site.append(i)
        elif freq == 0:
            remove_site.append(i)
        # multiallelic site
        elif len(counts) > 2:
            # This line is necessary, as an ancestral state can be the third most
            # common allele.
            del counts[site.ancestral_state]
            mutation_index = common_mutation_id(
                site, state=[counts.most_common(1)[0][0], site.ancestral_state]
            )
            keep_mutation.extend(mutation_index)
        # remove sites with 0 frequency
        else:
            keep_mutation.extend([mutation.id for mutation in site.mutations])

    tables = ts.dump_tables()
    drop_mutations(tables, keep_mutation)
    ts = tables.tree_sequence()
    down_sample_ts = ts.delete_sites(remove_site)
    return down_sample_ts


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

    seed = int(snakemake.params.individual_seed)

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
        pop_ts = maf_threshold(pop_ts, maf=0.05)

        with tempfile.TemporaryDirectory() as temp_dir:
            vcf_file_name = os.path.join(temp_dir, "analysis.vcf")

            with open(vcf_file_name, "w") as vcf_file:
                pop_ts.write_vcf(vcf_file)

            lddecay_command = [
                snakemake.params.poplddecay,
                "-InVCF",
                vcf_file_name,
                "-MAF",
                "0.05",
                "-OutStat",
                snakemake.output[f"{pop}_ld_decay"],
            ]

            subprocess.run(
                lddecay_command, check=True, stdout=sys.stderr, stderr=sys.stderr
            )


if __name__ == "__main__":
    main()
