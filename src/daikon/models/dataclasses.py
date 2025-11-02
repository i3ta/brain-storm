from dataclasses import dataclass
from typing import Self

import dill
import numpy as np
import numpy.typing as npt


class BaseDataClass:
    def save(self, filename: str):
        with open(filename, "wb") as f:
            dill.dump(self, f)

    @classmethod
    def load(cls, filename: str) -> Self:
        with open(filename, "rb") as f:
            return dill.load(f)


@dataclass
class TrainingResults(BaseDataClass):
    """
    Dataclass for storing data from training
    """

    model_name: str  # Model used
    loss: npt.NDArray[np.float32]  # Loss from each epoch of training

    # Training parameters
    batch_size: int
    learning_rate: int
    patience: int


@dataclass
class ModelEvaluation(BaseDataClass):
    """
    Dataclass for storing model evaluation results
    """

    model_name: str
    labels: list[str]
    y_pred: npt.NDArray[np.float32]
    y_true: npt.NDArray[np.float32]
    confusion: npt.NDArray[np.float32]
    f1: float
    fpr: npt.NDArray[np.float32]
    tpr: npt.NDArray[np.float32]
    auc_roc: float
