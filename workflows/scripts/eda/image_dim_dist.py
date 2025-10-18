import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from daikon.eda.dataclasses import DataAnalysisOutput


def plot_dims(type: str, input_file: str, ax: Axes):
    """
    Visualize the distribution of the image dimensions.
    """
    data = DataAnalysisOutput.load(input_file)
    dims = data.dims[:, :2]

    # Plot variables
    max_bar_count = 16

    # Get data for image dimensions
    dim_values, dim_counts = np.unique(dims, axis=0, return_counts=True)
    sort_idx = np.argsort(dim_counts)[::-1]
    dim_values = dim_values[sort_idx[:max_bar_count]]
    dim_counts = dim_counts[sort_idx[:max_bar_count]]
    bar_labels = [f"{w}x{h}" for w, h in dim_values]
    shown = min(max_bar_count, len(bar_labels))

    # Create plot
    bars = ax.bar(bar_labels, dim_counts)
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{int(height)}",
            ha="center",
            va="bottom",
        )

    ax.tick_params(axis="x", rotation=90)
    ax.set_ylabel("Frequency")
    ax.set_xlabel(f"Resolution (Width × Height)")
    ax.set_title(
        f"Image Resolution Distribution: {type.replace('_', ' ').title()}\n(n={dims.shape[0]:,}, top {shown} shown)"
    )


def main(
    input_files: list[str],
    output_file: str,
):
    """
    Visualize the distribution of the image dimensions.
    """
    types = ["glioma", "meningioma", "no_tumor", "pituitary"]

    # Create plot
    fig, axs = plt.subplots(2, 2, figsize=(16, 12))
    for i in range(len(input_files)):
        plot_dims(types[i], input_files[i], axs[i // 2, i % 2])
    fig.tight_layout()
    fig.savefig(output_file)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(
        input_files=[
            snakemake.input[0],
            snakemake.input[1],
            snakemake.input[2],
            snakemake.input[3],
        ],
        output_file=snakemake.output[0],
    )
