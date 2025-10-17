import os
import sys

import numpy as np
from tqdm import tqdm

from daikon.eda import DataAnalysisOutput, get_basic_metrics


def main(data: str, output: str, use_tqdm=False):
    """
    Perform basic exploratory data analysis on all images in one folder and save data as .pkl
    file

    Args:
        data (str): image data folder
        output (str): folder for output
    """
    image_extensions = (".png", ".jpg", ".jpeg")

    n = 0
    all_dims: list[tuple[int, int, int]] = []
    all_contr: list[float] = []

    for image in tqdm(os.listdir(data), disable=not use_tqdm):
        if image.lower().endswith(image_extensions):
            n += 1
            image_path = os.path.join(data, image)
            dims, contr = get_basic_metrics(image_path)
            all_dims.append(dims)
            all_contr.append(contr)

    data_output = DataAnalysisOutput(n, np.array(all_dims), np.array(all_contr))
    data_output.save(output)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(data=snakemake.input[0], output=snakemake.output[0])

elif __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            f"Expected 2 arguments, found {len(sys.argv) - 1}. Script should be called with data_folder, output_file"
        )
        exit(1)

    main(data=sys.argv[1], output=sys.argv[2], use_tqdm=True)
