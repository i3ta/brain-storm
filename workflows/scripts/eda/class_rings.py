import os

import numpy as np
import numpy.typing as npt

from daikon.eda.clustering import generate_circle_brightness
from daikon.eda.dataclasses import ImagePixelRings


def main(dir: str, output_file: str):
    """
    Generate histograms for all images in directory.

    Args:
        dir (str): image data folder
        output_file (str): output file path
    """
    image_extensions = (".png", ".jpg", ".jpeg")

    n = 0
    hist: dict[str, npt.NDArray[np.float64]] = {}
    for image in os.listdir(dir):
        if image.lower().endswith(image_extensions):
            n += 1
            image_path = os.path.join(dir, image)
            hist[image_path] = generate_circle_brightness(image_path, 32)

    data_output = ImagePixelRings(n, hist)
    data_output.save(output_file)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(dir=snakemake.input[0], output_file=snakemake.output[0])
