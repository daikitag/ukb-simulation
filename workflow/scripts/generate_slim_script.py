import string
import textwrap
import demes
from pathlib import Path

from utils import obtain_msprime_ratemap

from snakemake.script import snakemake as snk

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

    source("resources/demes-slim/demes.slim");
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
    sim.treeSeqOutput(trees_file);
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
    mutation_rate, recombination_map, demography_file,
    scaling_selection, seed
):
    recomb_rates, recomb_ends = msprime_rm_to_slim_rm(recombination_map)
    indent = 8 * " "
    recomb_rates_str = slim_array_string(recomb_rates, indent)
    recomb_ends_str = slim_array_string(recomb_ends, indent)

    demography = demes.load(demography_file)
    demography_file_path = Path(demography_file)
    json_demography_file = demography_file_path.with_suffix(".json")
    demes.dump(demography, filename=json_demography_file, format="json", simplified=False)

    final_slim_script =  string.Template(slim_script).substitute(
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
    chromosome = int(snk.params.chromosome)
    arm = snk.params.arm

    seed = int(snk.params.slim_seed)

    seed *= chromosome
    seed += 1 if arm == "p" else 0

    recombination_map, _ = obtain_msprime_ratemap(
        recombination_map_file=snk.input.recombination_map_file,
        position_file=snk.input.position_file,
        chromosome=chromosome+arm,
    )

    slim_script = slim_makescript(
        mutation_rate=float(snk.params.mu),
        recombination_map=recombination_map,
        demography_file=snk.input.demography,
        scaling_selection=float(snk.params.selection_scaling),
        seed=seed,
    )

    with open(snk.output[0], "w") as f:
        print(slim_script, file=f)

if __name__=="__main__":
    main()
