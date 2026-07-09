import sys

import demes
import msprime
import numpy as np
import pandas as pd
import pyslim
import tskit
from utils import obtain_msprime_ratemap


def subset_tree_seq(ts, selected_individuals):
    """Subset input tree sequence data based on the individual IDs.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Input tree sequence data.
    selected_individuals : numpy.array
        Numpy array of individual IDs that will be included in the output tree sequence
        data.

    Returns
    -------
    subset_ts : tskit.TreeSequence
        Subset of the input tree sequence data with the individuals specified in
        `selected_individuals`.

    Notes
    -----
    Since the simplification of the tree sequence data is conducted by using
    `filter_individuals=False`, the individual IDs will be preserved in the output
    tree sequence data. This allows the user to directly use the individual IDs to
    specify the individuals in the new output tree sequence data without converting
    the individual IDs.
    """
    selected_nodes = np.array([], dtype=int)
    for individual in selected_individuals:
        selected_nodes = np.concatenate(
            (selected_nodes, ts.individual(individual).nodes)
        )

    subset_ts = ts.simplify(selected_nodes, filter_individuals=False)

    return subset_ts


def convert_allele(ts):
    """Generate nucleotide alleles for input tree sequence data.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Input tree sequence data in SLiM alleles. This should be an output from SLiM
        simulation.

    Returns
    -------
    ts : tskit.TreeSequence
        Output tree sequence data in nucleotide alleles. The information of the input
        tree sequence data outside of alleles (e.g., individual and mutation
        information) are not modified.
    """
    ts = pyslim.generate_nucleotides(ts)
    ts = pyslim.convert_alleles(ts)

    return ts


def obtain_demes_demography(demes_demography):
    """Obtain msprime demography from demes YAML file.

    Parameters
    ----------
    demes_demography : str
        Path to demes YAML file.

    Returns
    -------
    demography : msprime.Demography
        A `Demography` instance of `msprime` that corresponds to the input demes
        YAML file.

    Notes
    -----
    The output of this function can be directly loaded into msprime's ancestry
    simulation as a demography input.
    """
    graph = demes.load(demes_demography)
    demography = msprime.Demography.from_demes(graph)

    return demography


def recapitate_shift(ts, demography, recombination_rate, left_position, random_seed):
    """Recapitate the input tree sequence data and shift its coordinates.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Input tree sequence file.
    demography: msprime.Demography
        Demography instance that will be used to conduct recapitation.
    recombination_rate : msprime.RateMap
        Recombination map that will be used to conduct recapitation.
    left_position : int
        Integer that will be used to shift the coordinate system.

    Returns
    -------
    ts : tskit.TreeSequence
        Output tree sequence data with recapitation and shifted coordinates.

    Notes
    -----
    Three steps are conducted in this function:

    (1) Recapitation
    It is required that the input tree sequence file has not coalesced completely
    (this should be the case when the tree sequence file is generated from SLiM
    simulation), and a neutral coalescent simulation is conducted to simulate
    ancestry information of individuals that are present in the input tree
    sequence data.

    (2) Simplification
    Simplification is conducted after recapitation to reduce nodes and mutations
    that are not present in the ancestral information of individuals that exist
    in present day. Note that this step is conducted after recapitation, as
    some of these internal nodes are necessary for recapitation.

    (3) Shift
    The coordinate of the output SLiM simulation is zero-based, so it is vital
    for us to shift the coordinate system of the output simulation.

    For recapitation and simplification, please refer to `pyslim` documentation,
    as it describes these concepts in detail.
    """
    ts = msprime.sim_ancestry(
        initial_state=ts,
        demography=demography,
        recombination_rate=recombination_rate,
        random_seed=random_seed,
    )

    ts = ts.simplify()

    ts = ts.shift(left_position)

    return ts


def main():
    """Recapitate and sample individuals.

    This function will first recapitate the input tree sequence data and shift
    its coordinate system to match the recombination map. Afterwards, it takes a
    subset of the tree sequence data based on the input individual ID dataframe.
    """
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    ts = tskit.load(snakemake.input.ts)

    chromosome = int(snakemake.params.chromosome)
    arm = snakemake.params.arm
    recapitate_seed = int(snakemake.params.recapitate_seed)

    # This is to make sure that we can use the same seed to conduct whole genome
    # simulation. Users do not need to specify a different seed for each chromosome.
    recapitate_seed *= chromosome
    recapitate_seed += 1 if arm == "p" else 0

    recombination_map, left_position = obtain_msprime_ratemap(
        recombination_map_file=snakemake.input.recombination_map_file,
        position_file=snakemake.input.position_file,
        chromosome=str(chromosome) + arm,
    )

    demography = obtain_demes_demography(snakemake.input.demography)

    ts = recapitate_shift(
        ts=ts,
        demography=demography,
        recombination_rate=recombination_map,
        left_position=left_position,
        random_seed=recapitate_seed,
    )

    individual_id_df = pd.read_csv(snakemake.input.individual_id)

    ts = subset_tree_seq(ts, individual_id_df["individual_id"].astype(int))

    ts = convert_allele(ts)

    ts.dump(snakemake.output.ts)


if __name__ == "__main__":
    main()
