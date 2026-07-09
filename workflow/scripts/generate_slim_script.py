import string
import sys
import textwrap
from pathlib import Path

import demes
from utils import obtain_msprime_ratemap

"""A brief explanation to the SLiM script below:

This SLiM script conducts a forward-time simulation where individuals are subject to
underdominant selection model. The selection coefficient of each mutation is randomly
selected from an input csv file with simulated values of s from Simons et al. (2025).

For this mutation model, a mutation has a selection coefficient of
-0.5 * `s` / `scaling_selection`, and a dominance coefficient of `scaling_selection`.
In SLiM, the relative fitness of an individual is modeled as
1, 1 + selection coefficient * dominance coefficient, and 1 + selection coefficient.
In an underdominant selection model, we are interested in modeling individual's fitness
as 1, 1 + selection_coefficient / 2, and 1. While it would be possible for us to
directly use this in SLiM simulation, the computational speed of SLiM simulation would
be the fastest when we simply model an individual's fitness by using the mutation's
selection coefficient and dominance coefficient. We can accomplish this selection
model by setting the `scaling_selection` parameter as a really small number.

The biggest limitation of running SLiM for a biobank-scale sample size is the required
RAM. The required RAM only gets extremely large at the final generations when we have
an explosive growth. To reduce the RAM, we implement two things: (1) setting the
mutation rate as 0 in the final generations, and (2) manually ask SLiM to run
simplication every five generations in the final 42 generations. Running simplification
in SLiM can result in reduced memory, but with a cost of increased computational time.
So we use SLiM's automatic simplification for the majority of generations and only run
repeated simplication at the final generations.

The output file name (tree_filename) is intentionally left as an undefined variable.
This will be defined in the snakemake pipeline.

Demes-slim is used to load the demography file into SLiM.

The below script is written as a template string so that we can input an arbitrary
number for each of these variables.
"""
slim_script = """
initialize() {

    initializeTreeSeq(timeUnit="generations");
    initializeMutationRate($mutation_rate);

    defineConstant("scaling_selection", $scaling_selection);
    df = readCSV("resources/simulated_s_001.csv", sep="\\r");
    defineConstant("simulated_s", df.getValue('simulated_s'));
    scriptForQTLs = "-0.5 * sample(simulated_s, 1) / scaling_selection;";
    initializeMutationType("m1", scaling_selection, "s", scriptForQTLs);

    defineConstant("recombination_rates", $recombination_rates);
    defineConstant("recombination_ends", $recombination_ends);
    initializeRecombinationRate(recombination_rates, recombination_ends);

    // g1 genomic element type: uses m1 for all mutations
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, $element_starts, $element_ends);

    setSeed($seed);

    source("demes-slim/demes.slim");
    defineConstant("SCALING_FACTOR", 1.0);
}

1 early() {
    model = demes_load("$demography_file", scaling_factor=SCALING_FACTOR, burn_in=0.0);
    demes_schedule_events(model);
    defineConstant("sim_end", model.getValue("end_time"));
    community.rescheduleScriptBlock(s1, sim_end, sim_end);
}

early() {
    if (community.tick == (sim_end - 40))
        sim.chromosome.setMutationRate(0);
}

late() {
    if ((community.tick >= (sim_end - 42)) &
        ((community.tick - (sim_end - 42)) % 5 == 0))
        sim.treeSeqSimplify();
}

s1 late() {
    sim.treeSeqSimplify();
    sim.treeSeqOutput(tree_filename);
}
"""


def slim_array_string(iterable, indent, width=80):
    """
    Format an array as a SLiM c() array and return as a line-wrapped string.
    """
    return (
        "c(\n"
        + textwrap.fill(
            ", ".join(map(str, iterable)),
            width=width,
            initial_indent=indent,
            subsequent_indent=indent,
        )
        + ")"
    )


def msprime_rm_to_slim_rm(recombination_map):
    """
    Taken from: https://github.com/popsim-consortium/stdpopsim/blob/8bc753eb9743531b4ba5205ecc75073a672711ff/stdpopsim/slim_engine.py#L884
    Convert recombination map from start position coords to end position coords.

    In SLiM, if ends[j-1] = a and ends[j] = b, then the recombination rate rates[j]
    applies to the links between a and b, i.e., to the links a:(a+1), (a+1):(a+2),
    ... (b-1):b. The tree sequence output by a SLiM simulation with L loci
    (i.e., positions 0, ..., L-1) will have sequence length equal to L, because
    intervals in tskit are open on the right, so the interval [0, L) does not
    include L.

    On the other hand, in msprime, a recombination rate map with some rate
    applied to the interval [x, y) will allow recombination events to the
    integers falling in [x, y); an event occuring at x will split x-1 from x,
    and so this implies recombination for the links from
    (x-1):x, x:(x+1), ..., (y-2):(y-1); this would correspond to ends of x-1
    and y-1 in SLiM.

    Note that this implies that the recombination rate that a msprime RateMap
    assigns to the interval [0, 1) has no effect in a discrete msprime
    simulation.
    """
    rates = recombination_map.rate.copy()
    # replace missing values with 0 recombination rate
    rates[recombination_map.missing] = 0
    ends = [int(pos) - 1 for pos in recombination_map.position]
    return rates, ends[1:]


def slim_makescript(
    mutation_rate,
    recombination_map,
    demography_file,
    scaling_selection,
    seed,
    slim_script=slim_script,
):
    """Generate SLiM script.

    This function uses the input template string of the SLiM script and substitutes
    the input varibles to generate the final SLiM script that can be used in SLiM
    simulation.

    Parameters
    ----------
    mutation_rate : float
        Mutation rate per base per generation.
    recombination_map : msprime.RateMap
        Recombination map in msprime's structure.
    demography_file : str
        File path to demes demography YAML file.
    selection_scaling : float
        Selection scaling parameter in SLiM simulation.
    seed : int
        Seed that will be used in SLiM simulation.
    slim_script : str
        A template string of SLiM script.

    Returns
    -------
    final_slim_script : str
        SLiM script that can be used in SLiM simulation.
    """
    recomb_rates, recomb_ends = msprime_rm_to_slim_rm(recombination_map)
    indent = 8 * " "
    recomb_rates_str = slim_array_string(recomb_rates, indent)
    recomb_ends_str = slim_array_string(recomb_ends, indent)

    demography = demes.load(demography_file)
    demography_file_path = Path(demography_file)
    json_demography_file = demography_file_path.with_suffix(".json")
    demes.dump(
        demography, filename=json_demography_file, format="json", simplified=False
    )

    final_slim_script = string.Template(slim_script).substitute(
        recombination_rates=recomb_rates_str,
        recombination_ends=recomb_ends_str,
        element_starts=0,
        element_ends=recomb_ends[-1],
        mutation_rate=mutation_rate,
        demography_file=json_demography_file,
        scaling_selection=scaling_selection,
        seed=seed,
    )

    return final_slim_script


def main():
    """Generates SLiM script.

    This function generates a SLiM script that can be used in conducting SLiM
    simulation. A template string is used to describe the initial SLiM script, as we
    would like to use custom Python functions to input the recombination map in a
    correct format to SLiM. The final SLiM script can be directly used in SLiM, so it
    is saved for later usage.

    The function to load the recombination map into msprime and convert it to SLiM
    format is adapted from the simulation pipeline in `stdpopsim`.
    """
    sys.stderr = open(snakemake.log[0], "w", buffering=1)
    chromosome = int(snakemake.params.chromosome)
    arm = snakemake.params.arm

    seed = int(snakemake.params.slim_seed)

    # This is used to set the seed for each chromosome and arm as a different
    # interger
    seed *= chromosome
    seed += 1 if arm == "p" else 0

    recombination_map, _ = obtain_msprime_ratemap(
        recombination_map_file=snakemake.input.recombination_map_file,
        position_file=snakemake.input.position_file,
        chromosome=str(chromosome) + arm,
    )

    slim_script = slim_makescript(
        mutation_rate=float(snakemake.params.mu),
        recombination_map=recombination_map,
        demography_file=snakemake.input.demography,
        scaling_selection=float(snakemake.params.selection_scaling),
        seed=seed,
    )

    with open(snakemake.output.slim_script, "w") as f:
        print(slim_script, file=f)


if __name__ == "__main__":
    main()
