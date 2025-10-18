import matplotlib.pyplot as plt
import numpy as np

from daikon.eda.dataclasses import DataAnalysisOutput


def main(
    glioma: str,
    meningioma: str,
    no_tumor: str,
    pituitary: str,
    output_img: str,
):
    """
    Collect the data for all types of folders and output a figure summary of
    data.
    """
    data = {
        "glioma": DataAnalysisOutput.load(glioma),
        "meningioma": DataAnalysisOutput.load(meningioma),
        "no_tumor": DataAnalysisOutput.load(no_tumor),
        "pituitary": DataAnalysisOutput.load(pituitary),
    }

    # Plot variables
    bar_width = 0.2
    bar_group_count = 8
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, 4))

    # Get data for image dimensions
    all_dims = np.vstack([d.dims for d in data.values()])
    dim_values, dim_counts = np.unique(all_dims[:, :2], axis=0, return_counts=True)
    sort_idx = np.argsort(dim_counts)[::-1]
    dim_values = dim_values[sort_idx[:bar_group_count]]
    bar_labels = [f"{w}x{h}" for w, h in dim_values]

    # Get data for contrast
    all_contr = [d.contr for d in data.values()]
    contr_labels = data.keys()

    # Create plot
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(12, 12))
    n = sum([d.n for d in data.values()])

    # Plot image dimensions
    x = np.arange(len(bar_labels))
    for i, (category, data) in enumerate(data.items()):
        dims = data.dims[:, :2]
        dim_counts = []
        for dv in dim_values:
            count = np.sum(np.all(dims == dv, axis=1))
            dim_counts.append(count)

        offset = bar_width * i
        color = colors[i % len(colors)]
        rects = ax1.bar(
            x + offset,
            dim_counts,
            bar_width,
            label=f"{category} (n = {dims.shape[0]})",
            color=color,
        )
        ax1.bar_label(rects, padding=3)

    ax1.set_title(f"Common Image Dimensions (n = {n})")
    ax1.set_ylabel("Number of Images")
    ax1.set_xlabel("Image Dimensions (px)")
    ax1.set_xticks(x + bar_width, bar_labels)
    ax1.legend(loc="upper right")

    # Plot violin plots
    ax2.violinplot(all_contr, showextrema=False)
    ax2.boxplot(
        all_contr,
        widths=0.1,
        patch_artist=True,
        boxprops=dict(facecolor="white", color="black"),
        medianprops=dict(color="red", linewidth=1.5),
        whiskerprops=dict(color="black"),
        capprops=dict(color="black"),
    )
    for i, contr in enumerate(all_contr):
        sample_size = max(1, int(len(contr) * 0.05))
        sample_indices = np.random.choice(len(contr), size=sample_size, replace=False)
        sampled_contr = contr[sample_indices]

        x_jitter = np.random.normal(i + 1, 0.03, size=len(sampled_contr))
        ax2.scatter(x_jitter, sampled_contr, color="black", s=8, alpha=0.6, zorder=3)

    ax2.set_xticks(np.arange(1, len(contr_labels) + 1), labels=contr_labels)
    ax2.set_xlim(0.2, len(contr_labels) + 0.8)
    ax2.set_title(f"Normalized Grayscale Image Contrasts (n = {n})")
    ax2.set_xlabel("Category")
    ax2.set_ylabel("Normalized Contrast")

    # Save figure
    fig.tight_layout()
    fig.savefig(output_img)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(
        glioma=snakemake.input[0],
        meningioma=snakemake.input[1],
        no_tumor=snakemake.input[2],
        pituitary=snakemake.input[3],
        output_img=snakemake.output[0],
    )
