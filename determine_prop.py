import sys

import numpy as np
import pandas as pd
import yaml


def main():
    l_inferred = float(sys.argv[1])
    recom_map = pd.read_csv("resources/recombination_map/recombination_position.csv")
    genome_length = np.sum(recom_map["right"] - recom_map["left"])

    with open("config/config.yaml") as yml:
        config = yaml.safe_load(yml)

    causal_mu = float(config["mu"])

    prop = (l_inferred / genome_length) * (1.25e-8 / causal_mu)

    print("Prop is:", prop)


if __name__ == "__main__":
    main()
