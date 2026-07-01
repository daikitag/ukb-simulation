import pandas as pd
import msprime


def obtain_msprime_ratemap(recombination_map_file, position_file, chromosome):
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
