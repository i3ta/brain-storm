import matplotlib.pyplot as plt
import numpy as np

from daikon.eda.dataclasses import DataAnalysisOutput


def main(type: str, input_file: str, output_file: str):
    """
    Visualize the distribution of the image dimensions.
    """
    data = DataAnalysisOutput.load(input_file)
    dims = data.dims[:, :2]

    # Plot variables
    bar_count = 32

    # Get data for image dimensions
    dim_values, dim_counts = np.unique(dims, axis=0, return_counts=True)
    sort_idx = np.argsort(dim_counts)[::-1]
    dim_values = dim_values[sort_idx[:bar_count]]
    dim_counts = dim_counts[sort_idx[:bar_count]]
    bar_labels = [f"{w}x{h}" for w, h in dim_values]

    # Create plot
    plt.figure(figsize=(10, 6))
    plt.bar(bar_labels, dim_counts)

    plt.xticks(rotation=90)
    plt.ylabel("Frequency")
    plt.xlabel(f"Resolution (Width × Height)")
    plt.title(
        f"Image Resolution Distribution: {type.replace('_', ' ').title()}\n(n={dims.shape[0]:,}, top {bar_count} shown)"
    )
    plt.tight_layout()
    plt.savefig(output_file)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(
        type=snakemake.wildcards.type,
        input_file=snakemake.input[0],
        output_file=snakemake.output[0],
    )
