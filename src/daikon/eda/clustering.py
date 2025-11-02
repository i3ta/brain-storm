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
        hist[br] = np.count_nonzero(image_flat == br) / image.shape[1]

    return hist


def generate_circle_brightness(image_path: str, rings: int) -> npt.NDArray[np.float64]:
    """
    Generate the circle pixel brightness for a given image.

    The image is divided into rings, and the average pixel brightness is
    calculated for each ring. These brightnesses are used in our clustering
    """
    image = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    if image is None:
        raise RuntimeError("Error reading image")

    center = (image.shape[0] / 2, image.shape[1] / 2)
    x_coords, y_coords = np.indices((image.shape[0], image.shape[1]))
    dist = np.sqrt((x_coords - center[0]) ** 2 + (y_coords - center[0]) ** 2)
    width = max(image.shape[0], image.shape[1]) / rings
    max_br = np.max(image)

    br = np.zeros(rings)
    for r in range(rings):
        r_min = r * width
        r_max = (r + 1) * width

        ring_mask = (dist >= r_min) & (dist < r_max)
        pixels = image[ring_mask]

        if len(pixels) > 0:
            br[r] = np.mean(pixels, dtype=np.float64) / max_br

    return br
