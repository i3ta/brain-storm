from typing import Any

import cv2
import numpy as np


def generate_seed(seed: int, seed_str: str):
    """
    Generate a seed from a seed number and string
    """
    return sum(ord(c) for c in seed_str) + seed


def duplicate_image_augment(source: str, target: str, seed: int, cfg: dict[str, Any]):
    """
    Duplicate an image and augment it to produce a new image
    """
    rng = np.random.default_rng(generate_seed(seed, target))
    img = cv2.imread(source, cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise RuntimeError("Error reading image")

    rows, cols = img.shape

    # Augmentation parameters
    rotation = rng.uniform(cfg["rotation_low"], cfg["rotation_high"])
    brightness_mult = rng.uniform(cfg["brightness_low"], cfg["brightness_high"])
    brightness_shift = rng.uniform(
        cfg["brightness_shift_low"], cfg["brightness_shift_high"]
    )
    flip = rng.uniform()
    noise = rng.normal() * np.sqrt(cfg["noise_variance"])

    # Apply augmentation
    M = cv2.getRotationMatrix2D((cols / 2, rows / 2), rotation, 1)
    img = cv2.warpAffine(img, M, (cols, rows))
    img = cv2.convertScaleAbs(img, alpha=brightness_mult, beta=brightness_shift)
    if flip > 0.5:
        img = cv2.flip(img, 1)
    img = np.clip(img + noise, 0, 255).astype(np.uint8)

    # Save image
    cv2.imwrite(target, img)


def resize_normalize(image_path: str, dest_path: str, w: int, h: int):
    """
    Resize the image to given dimensions
    """
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise RuntimeError("Error reading image")

    orig_h, orig_w = image.shape[:2]

    # Compute scale so that the smaller side matches the desired size
    scale = max(w / orig_w, h / orig_h)

    # Resize with the computed scale
    new_w = max(int(orig_w * scale), w)
    new_h = max(int(orig_h * scale), h)
    resized = cv2.resize(image, (new_w, new_h), interpolation=cv2.INTER_AREA)

    # Compute the crop coordinates (center crop)
    start_x = (new_w - w) // 2
    start_y = (new_h - h) // 2

    cropped = resized[start_y : start_y + h, start_x : start_x + w]

    # Normalize to 0–255 based on min/max of the cropped image
    min_val = cropped.min()
    max_val = cropped.max()

    if max_val > min_val:  # avoid divide-by-zero
        normalized = (cropped - min_val) * (255.0 / (max_val - min_val))
    else:
        normalized = np.zeros_like(cropped, dtype=np.float32)

    cv2.imwrite(dest_path, normalized)
