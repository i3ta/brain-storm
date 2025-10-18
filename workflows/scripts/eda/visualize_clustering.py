import os

import cv2
import matplotlib.pyplot as plt
import numpy as np

from daikon.eda.clustering import generate_histogram
from daikon.eda.dataclasses import ImagePixelClustering


def main(
    input_clusters: list[str],
    data_dirs: str,
    k: int,
    output_line: str,
    output_examples: str,
):
    """
    Visualize the clustering at a given number of clusters.

    Args:
        input_clusters (list[str]): List of paths for all clusters
        data (str): Folder for data
        k (int): Target number of clusters
        output_file (str): Output file name
    """
    clusters = np.arange(1, 19)
    distortion = np.zeros(len(input_clusters))
    target_kmeans = None
    for i in range(len(input_clusters)):
        data = ImagePixelClustering.load(input_clusters[i])
        distortion[i] = data.kmeans.inertia_
        if data.k == k:
            target_kmeans = data.kmeans

    if target_kmeans == None:
        raise RuntimeError("Target number of clusters not found")

    # Set up line plot
    fig1 = plt.figure(figsize=(10, 8))
    ax1 = fig1.add_subplot(111)
    ax1.plot(clusters, distortion)
    ax1.axvline(x=k, color="r", linestyle="--", label=f"k={k}")

    k_index = k - 1
    distortion_at_k = distortion[k_index]
    y_range = ax1.get_ylim()[1] - ax1.get_ylim()[0]
    vertical_offset = 0.05 * y_range

    ax1.text(
        k,
        distortion_at_k + vertical_offset,
        f"  Distortion: {distortion_at_k:.2e}",
        verticalalignment="center",
        color="r",
    )

    ax1.set_ylabel("Distortion")
    ax1.set_xlabel("Clusters")
    ax1.legend()
    fig1.tight_layout()
    fig1.savefig(output_line)

    # Set up example images plot
    row_labels = ["glioma", "meningioma", "no_tumor", "pituitary"]
    fig2, axs = plt.subplots(4, k, figsize=(k * 2, 8))

    # Select data for examples
    image_extensions = (".png", ".jpg", ".jpeg")
    found = set()  # Classes that we already have images for
    for i, dir in enumerate(data_dirs):
        added = 0
        for image in os.listdir(dir):
            if image.lower().endswith(image_extensions):
                img_path = os.path.join(dir, image)
                hist = generate_histogram(img_path)
                res = target_kmeans.predict([hist])
                if (i, res[0]) not in found:
                    found.add((i, res[0]))
                    img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
                    axs[i, res[0]].imshow(img, cmap="gray")
                    added += 1
                if added == k:
                    break

    for i, label in enumerate(row_labels):
        axs[i, 0].set_ylabel(label)
    for j, label in enumerate(range(k)):
        axs[0, j].set_title(label)
    for ax in axs.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["left"].set_visible(False)

    fig2.tight_layout()
    fig2.savefig(output_examples)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    n = len(snakemake.input)

    main(
        input_clusters=snakemake.input[:-4],
        data_dirs=snakemake.input[-4:],
        k=int(snakemake.wildcards.k),
        output_line=snakemake.output[0],
        output_examples=snakemake.output[1],
    )
