import sys

import matplotlib.pyplot as plt
import pandas as pd


def plot_pca(yri_df, ceu_df, jpt_df, chb_df, file_output, chromosome, arm):
    # point size, alpha values, axis font size
    size = 15
    alp = 0.3
    axis_fsize = 15

    fig = plt.figure(figsize=(15, 4))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    groups = {
        "YRI": (yri_df, "C1"),
        "CEU": (ceu_df, "C0"),
        "JPT": (jpt_df, "C2"),
        "CHB": (chb_df, "C3"),
    }

    for label, (df, color) in groups.items():
        ax1.scatter(
            df["PC1"],
            df["PC2"],
            marker=".",
            s=size,
            alpha=alp,
            label=label,
            color=color,
            rasterized=True,
        )

    # Hide the right and top spines
    ax1.spines["right"].set_visible(False)
    ax1.spines["top"].set_visible(False)

    legend = ax1.legend(fontsize=10, markerscale=3)
    for h in legend.legend_handles:
        h.set_alpha(1)

    ax1.set_xlabel("PC1", fontsize=axis_fsize)
    ax1.set_ylabel("PC2", fontsize=axis_fsize)

    for label, (df, color) in groups.items():
        ax2.scatter(
            df["PC3"],
            df["PC4"],
            marker=".",
            s=size,
            alpha=alp,
            label=label,
            color=color,
            rasterized=True,
        )

    # Hide the right and top spines
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)

    legend = ax2.legend(fontsize=10, markerscale=3)
    for h in legend.legend_handles:
        h.set_alpha(1)

    ax2.set_xlabel("PC3", fontsize=axis_fsize)
    ax2.set_ylabel("PC4", fontsize=axis_fsize)

    for label, (df, color) in groups.items():
        ax3.scatter(
            df["PC5"],
            df["PC6"],
            marker=".",
            s=size,
            alpha=alp,
            label=label,
            color=color,
            rasterized=True,
        )

    # Hide the right and top spines
    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)

    legend = ax3.legend(fontsize=10, markerscale=3)
    for h in legend.legend_handles:
        h.set_alpha(1)

    ax3.set_xlabel("PC5", fontsize=axis_fsize)
    ax3.set_ylabel("PC6", fontsize=axis_fsize)

    plt.suptitle(f"Chromosome {chromosome}{arm}: PCA", fontsize=axis_fsize * 1.2)

    fig.savefig(file_output, bbox_inches="tight", format="pdf", dpi=600)


def main():
    sys.stderr = open(snakemake.log[0], "w", buffering=1)

    pcs_df = pd.read_csv(snakemake.input.pcs, delimiter="\t")

    yri_df = pcs_df[pcs_df["IID"].str.contains("YRI")]
    ceu_df = pcs_df[pcs_df["IID"].str.contains("CEU")]
    jpt_df = pcs_df[pcs_df["IID"].str.contains("JPT")]
    chb_df = pcs_df[pcs_df["IID"].str.contains("CHB")]

    plot_pca(
        yri_df=yri_df,
        ceu_df=ceu_df,
        jpt_df=jpt_df,
        chb_df=chb_df,
        file_output=snakemake.output,
        chromosome=snakemake.params.chromosome,
        arm=snakemake.params.arm,
    )


if __name__ == "__main__":
    main()
