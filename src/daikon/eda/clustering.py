import cv2
import numpy as np
import numpy.typing as npt


def generate_histogram(image_path: str) -> npt.NDArray[np.int32]:
    """
    Generate the pixel brightness histogram for a given image.
    """
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise RuntimeError("Error reading image")

    image_flat = image.reshape(1, -1)
    hist = np.zeros(256, dtype=np.int32)
    for br in range(256):
        hist[br] = np.count_nonzero(image_flat == br)

    return hist
