import msprime
import pandas as pd


def obtain_msprime_ratemap(recombination_map_file, position_file, chromosome):
    """Obtain msprime rate map.

    This function loads the recombination map file in HapMap format and returns a
    RateMap object in msprime that can be used in msprime ancestry simulation.

    The `position_file` input is used to remove the centromeres from the simulation.
    The left position is given as an output as well, because it is used to preserve
    the coordinate system in the final simulation output.
    """
    position_df = pd.read_csv(position_file)
    recombination_map = msprime.RateMap.read_hapmap(
        recombination_map_file,
        rate_col=2,
    )
    recom_position = position_df[position_df.chromosome == chromosome]
    recombination_map = recombination_map.slice(
        left=recom_position.left.item(),
        right=recom_position.right.item(),
        trim=True,
    )

    return recombination_map, recom_position.left.item()
