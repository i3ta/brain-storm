import os

import numpy as np
import numpy.typing as npt

from daikon.eda.clustering import generate_histogram
from daikon.eda.dataclasses import ImagePixelHistogram


def main(dir: str, output_file: str):
    """
    Generate histograms for all images in directory.

    Args:
        dir (str): image data folder
        output_file (str): output file path
    """
    image_extensions = (".png", ".jpg", ".jpeg")

    n = 0
    hist: dict[str, npt.NDArray[np.int32]] = {}
    for image in os.listdir(dir):
        if image.lower().endswith(image_extensions):
            n += 1
            image_path = os.path.join(dir, image)
            hist[image_path] = generate_histogram(image_path)

    data_output = ImagePixelHistogram(n, hist)
    data_output.save(output_file)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(dir=snakemake.input[0], output_file=snakemake.output[0])
