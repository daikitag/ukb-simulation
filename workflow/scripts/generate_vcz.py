import tskit
import pandas as pd
import msprime
import numpy as np
import collections
import bio2zarr.tskit as ts2z
import tszip

from snakemake.script import snakemake as snk

def drop_mutations(tables, indexes_of_mutations_to_keep):
    m = len(tables.mutations)
    tables.mutations.parent = np.zeros(m, dtype=np.int32) - 1 # null the parent column
    select = np.zeros(m, dtype=np.bool)
    select[indexes_of_mutations_to_keep] = True
    tables.mutations.keep_rows(select)
    tables.compute_mutation_parents()

def common_mutation_id(site, state):
    mutation_list = []
    for m in site.mutations:
        if m.derived_state in state:
            mutation_list.append(m.id)
    return mutation_list

def _get_next_id(ts):
    max_id = -1
    for mut in ts.mutations():
        for d in mut.derived_state.split(","):
            max_id = max(max_id, int(d))
    return max_id + 1

def count_site_alleles(ts, tree, site):
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
            del counts[site.ancestral_state]
            mutation_index = common_mutation_id(site, state=[counts.most_common(1)[0][0],site.ancestral_state])
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
        selected_nodes = np.concatenate((selected_nodes, ts.individual(individual).nodes))

    return ts.simplify(selected_nodes)

def main():
    ts = tszip.load(snk.input.ts)

    ts = maf_threshold(ts, maf=float(snk.params.maf))

    individual_id_df = pd.read_csv(snk.input.individual_id)
    
    model_mapping = ts.map_to_vcf_model(
        individuals=individual_id_df["individual_id"],
        individual_names=individual_id_df["plink_id"],
        contig_id=str(snk.params.chromosome),
    )

    ts2z.convert(
        ts,
        vcz_path=snk.output.vcz,
        worker_processes=int(snk.threads),
        model_mapping=model_mapping,
    )


if __name__ == '__main__':
    main()