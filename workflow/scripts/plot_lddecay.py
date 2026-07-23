import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def ld_decay_meanbin(df, bin1=10, bin2=100, break_bp=100):
    short_mask = df["#Dist"] < break_bp
    bin_size = np.where(short_mask, bin1, bin2)

    df["bin_size"] = bin_size
    df["bin_idx"] = np.floor((df["#Dist"] - 0.1) / df["bin_size"]).astype(int)

    df["x_bp"] = (df["bin_idx"] + 1) * df["bin_size"]
    df["r2_x_count"] = df["Mean_r^2"] * df["NumberPairs"]

    binned = (
        df.groupby(["bin_size", "bin_idx", "x_bp"], as_index=False)
        .agg(
            sum_r2=("r2_x_count", "sum"),
            n_pairs=("NumberPairs", "sum"),
        )
        .sort_values("x_bp")
    )

    binned["mean_r2"] = binned["sum_r2"] / binned["n_pairs"]
    return binned[["x_bp", "mean_r2", "n_pairs"]]


def plot_lddecay(sim_df, thousand_df, chromosome, file_name):
    i = 1
    fig = plt.figure(figsize=(14, 10))
    for pop in ["ceu", "yri", "chb", "jpt"]:
        plt.subplot(2, 2, i)
        plt.plot(
            thousand_df[pop]["x_bp"] / 1000,
            thousand_df[pop]["mean_r2"],
            lw=2,
            label="1000 Genomes Project",
        )
        plt.plot(
            sim_df[pop]["x_bp"] / 1000,
            sim_df[pop]["mean_r2"],
            lw=2,
            label="Simulation",
        )
        plt.legend()
        plt.title(f"Population: {pop}")
        plt.xlabel("Distance (Kb)")
        plt.ylabel("Mean r2")
        i += 1

    plt.subplots_adjust(hspace=0.3, wspace=0.2)
    fig.suptitle(f"LD Decay Plot: Chromosome {chromosome}", fontsize=20)

    plt.savefig(file_name, bbox_inches="tight")


def main():
    thousand_df = {}
    sim_df = {}

    for pop in ["ceu", "yri", "chb", "jpt"]:
        sim_original_df = pd.read_csv(snakemake.input[f"{pop}_sim_ld_decay"], sep="\t")
        thousand_original_df = pd.read_csv(
            snakemake.input[f"{pop}_thousand_ld_decay"], sep="\t"
        )

        thousand_df[pop] = ld_decay_meanbin(sim_original_df)
        sim_df[pop] = ld_decay_meanbin(thousand_original_df)

    chromosome = str(snakemake.params.chromosome)
    arm = snakemake.params.arm

    plot_lddecay(
        sim_df=sim_df,
        thousand_df=thousand_df,
        chromosome=chromosome + arm,
        file_name=snakemake.output.plot,
    )


if __name__ == "__main__":
    main()
