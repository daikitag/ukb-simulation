import collections
import sys

import bio2zarr.tskit as ts2z
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


def main():
    """Generates VCZ.

    This function first subsets the input tree sequence data based on MAF and processes
    multi-allelic sites, such that all sites in the tree sequence data are biallelic.
    This is to make sure that no errors are raised when generating PLINK files from
    VCZ.

    The `plink_id` column in individual ID dataframe is used as individual names in
    VCZ, and it will be also be used in PLINK files.
    """
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    ts = tszip.load(snakemake.input.ts)

    ts = maf_threshold(ts, maf=float(snakemake.params.maf))

    individual_id_df = pd.read_csv(snakemake.input.individual_id)

    model_mapping = ts.map_to_vcf_model(
        individuals=individual_id_df["individual_id"],
        individual_names=individual_id_df["plink_id"],
        contig_id=str(snakemake.params.chromosome),
    )

    ts2z.convert(
        ts,
        vcz_path=snakemake.output.vcz,
        worker_processes=int(snakemake.threads),
        model_mapping=model_mapping,
    )


if __name__ == "__main__":
    main()
