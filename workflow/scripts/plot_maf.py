import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    sim_data = {}
    thousand_data = {}

    chromosome = str(snakemake.params.chromosome)
    arm = snakemake.params.arm

    for pop in ["ceu", "yri", "chb", "jpt"]:
        sim_df = pd.read_csv(snakemake.input[f"{pop}_sim_maf"])
        thousand_df = pd.read_csv(
            snakemake.input[f"{pop}_thousand_maf"], header=None, names=["MAF"], sep="\t"
        )

        sim_df["plot_MAF"] = sim_df["MAF"].apply(lambda x: x if x < 0.5 else 1 - x)
        thousand_df["plot_MAF"] = thousand_df["MAF"].apply(
            lambda x: x if x < 0.5 else 1 - x
        )

        sim_data[pop] = sim_df.plot_MAF
        thousand_data[pop] = thousand_df.plot_MAF

    boundary = [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
    groups = ["10-15%", "15-20%", "20-25%", "25-30%", "30-40%", "40-50%"]

    fig = plt.figure(figsize=(12, 10))
    i = 1
    for pop in ["ceu", "yri", "chb", "jpt"]:
        ax = fig.add_subplot(2, 2, i)
        sim_plot, _ = np.histogram(sim_data[pop], bins=boundary)
        thousand_plot, _ = np.histogram(thousand_data[pop], bins=boundary)

        afs_plot = pd.DataFrame(
            {
                "Simulation": sim_plot,
                "1000 Genomes Project": thousand_plot,
            },
            index=groups,
        )

        afs_plot.plot(
            kind="bar",
            stacked=False,
            xlabel="Minor Allele Frequency",
            ylabel="Frequency of SNP",
            ax=ax,
            title=f"Population: {pop}, {np.sum(thousand_plot) / np.sum(sim_plot)}",
        )
        i += 1
    fig.suptitle(
        f"Chromosome {chromosome + arm}: Allele Frequency Spectrum", y=0.95, fontsize=15
    )
    plt.subplots_adjust(wspace=0.3, hspace=0.4)
    plt.savefig(snakemake.output.plot, bbox_inches="tight")


if __name__ == "__main__":
    main()
