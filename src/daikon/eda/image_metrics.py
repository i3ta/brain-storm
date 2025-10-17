import cv2
import numpy as np


def get_basic_metrics(img: str):
    """
    Get basic metadata about an image in the dataset.

    Args:
        img (str): Filename of the image to analyze

    Returns:
        (w, h, ch), contr: image dimensions and number of channels, and
        contrast
    """
    image = cv2.imread(img)
    dims = image.shape
    image_float = image.astype(np.float32) / 255.0
    grayscale = cv2.cvtColor(image_float, cv2.COLOR_BGR2GRAY)
    contr = grayscale.std()

    return dims, contr
