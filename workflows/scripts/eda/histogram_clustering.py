import numpy as np
from sklearn.cluster import KMeans

from daikon.eda.dataclasses import ImagePixelClustering, ImagePixelHistogram


def main(input_files: list[str], k: int, output_file: str):
    """
    Perform clustering on k clusters.

    Args:
        input_files (list[str]): List of input histograms to cluster on.
        k (int): Number of clusters
        output_file (str): Output file path
    """
    histograms = []
    for file in input_files:
        data: ImagePixelHistogram = ImagePixelHistogram.load(file)
        histograms.extend(data.hist.values())

    X = np.vstack(histograms)
    kmeans = KMeans(n_clusters=k).fit(X)
    output = ImagePixelClustering(k, kmeans)
    output.save(output_file)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(
        input_files=snakemake.input,
        k=int(snakemake.wildcards.k),
        output_file=snakemake.output[0],
    )
