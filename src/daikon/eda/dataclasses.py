from dataclasses import dataclass
from typing import Self

import dill
import numpy as np
import numpy.typing as npt
from sklearn.cluster import KMeans


class BaseDataClass:
    def save(self, filename: str):
        with open(filename, "wb") as f:
            dill.dump(self, f)

    @classmethod
    def load(cls, filename: str) -> Self:
        with open(filename, "rb") as f:
            return dill.load(f)


@dataclass
class DataAnalysisOutput(BaseDataClass):
    """
    Dataclass for outputs of exporatory data analysis.
    """

    n: int  # Number of images total
    dims: npt.NDArray[np.int32]  # Dimensions of each image (width, height, channels)
    contr: npt.NDArray[np.float32]  # Root mean square contrast of each image


@dataclass
class ImagePixelHistogram(BaseDataClass):
    """
    Dataclass for image pixel brightness histogram.
    """

    n: int  # Number of images total
    hist: dict[
        str, npt.NDArray[np.int32]
    ]  # Dictionary of image name to pixel brightness histogram


@dataclass
class ImagePixelClustering(BaseDataClass):
    """
    Dataclass for clustered image pixel brightness
    """

    k: int  # Number of clusters
    kmeans: KMeans  # KMeans result
