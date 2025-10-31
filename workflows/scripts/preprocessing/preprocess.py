import os
import shutil
from typing import Any

import numpy as np

from daikon.preprocessing import duplicate_image_augment, resize_normalize


def main(src: str, dest: str, cfg: dict[str, Any], seed: int):
    """
    Preprocess images.

    Args:
        src (str): Source directory of images
        dest (str): Destination directory of images
        cfg (dict[str, Any]): Snakemake config file
    """
    rng = np.random.default_rng(seed)
    image_extensions = (".png", ".jpg", ".jpeg")
    image_dims = (cfg["image_width"], cfg["image_height"])

    # Folders
    dest_tmp = os.path.join(dest, "tmp")
    dest_test = os.path.join(dest, "test")
    dest_train = os.path.join(dest, "train")

    os.makedirs(dest_tmp, exist_ok=True)
    os.makedirs(dest_test, exist_ok=True)
    os.makedirs(dest_train, exist_ok=True)

    # Resize images and get total number of images per class
    dirs = os.listdir(src)
    count = {}
    classes = []
    for dir in dirs:
        src_dir = os.path.join(src, dir)
        dest_dir = os.path.join(dest_tmp, dir)
        os.makedirs(dest_dir, exist_ok=True)
        if not os.path.isdir(src_dir):
            continue

        img_count = 0
        classes.append(dir)
        for img in os.listdir(src_dir):
            src_path = os.path.join(src_dir, img)
            if img.lower().endswith(image_extensions):
                # Is an image
                img_count += 1
                _, ext = os.path.splitext(img)

                dest_path = os.path.join(dest_dir, f"{dir}_{img_count:04d}{ext}")
                resize_normalize(src_path, dest_path, image_dims[0], image_dims[1])

        count[dir] = img_count

    if "image_count" in cfg:
        target_count = cfg["image_count"]
    else:
        target_count = max(count.values())  # Choose maximum class as target count

    # Augment images to create more data
    for c in classes:
        # Get all images for class
        src_dir = os.path.join(dest_tmp, c)
        images = [
            img for img in os.listdir(src_dir) if img.lower().endswith(image_extensions)
        ]

        images = rng.choice(images, size=(target_count - count[c]), replace=True)
        for img in images:
            count[c] += 1
            _, ext = os.path.splitext(img)

            src_path = os.path.join(src_dir, img)
            dest_path = os.path.join(src_dir, f"{c}_{count[c]:04d}{ext}")

            duplicate_image_augment(src_path, dest_path, seed, cfg)

    # Separate data into train and test
    test_count = int(target_count * cfg["test_prop"])
    for c in classes:
        # Get all images for class
        src_dir = os.path.join(dest_tmp, c)
        images = np.array(
            [
                img
                for img in os.listdir(src_dir)
                if img.lower().endswith(image_extensions)
            ]
        )

        # Get indices for testing and training
        test_idx = rng.choice(
            np.arange(0, target_count), size=test_count, replace=False
        )
        test_img = images[test_idx]

        # Move files into test and train folders
        dest_dir = os.path.join(dest_test, c)
        os.makedirs(dest_dir, exist_ok=True)
        for img in test_img:
            src_path = os.path.join(src_dir, img)
            shutil.move(src_path, dest_dir)

        dest_dir = os.path.join(dest_train, c)
        os.makedirs(dest_dir, exist_ok=True)
        for img in os.listdir(src_dir):
            src_path = os.path.join(src_dir, img)
            shutil.move(src_path, dest_dir)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(
        src=snakemake.input[0],
        dest=snakemake.output[0],
        cfg=snakemake.config["preprocessing"],
        seed=int(snakemake.config["seed"]),
    )
