#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main() -> None:
    input_path = Path("results/chr_feature_counts.tsv")
    output_path = Path("results/chr_feature_counts_plot.png")

    df = pd.read_csv(input_path, sep="\t")
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)

    # Keep the table order (already sorted by gene_per_Mb descending).
    chrom_order = df["chrom"].tolist()
    melted = df.melt(
        id_vars=["chrom"],
        value_vars=["n_gene", "n_exon_unique", "n_tRNA", "n_snoRNA"],
        var_name="feature",
        value_name="count",
    )
    melted["chrom"] = pd.Categorical(melted["chrom"], categories=chrom_order, ordered=True)

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), constrained_layout=True)

    sns.barplot(
        data=df,
        x="chrom",
        y="gene_per_Mb",
        color="#2a9d8f",
        ax=axes[0],
    )
    axes[0].set_title("Gene Density by Chromosome")
    axes[0].set_xlabel("")
    axes[0].set_ylabel("Genes per Mb")
    axes[0].tick_params(axis="x", rotation=45)

    sns.barplot(
        data=melted,
        x="chrom",
        y="count",
        hue="feature",
        ax=axes[1],
    )
    axes[1].set_title("Feature Counts by Chromosome")
    axes[1].set_xlabel("Chromosome")
    axes[1].set_ylabel("Count")
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].legend(title="Feature")

    fig.savefig(output_path, dpi=300)
    print(f"Saved plot to: {output_path}")


if __name__ == "__main__":
    main()
