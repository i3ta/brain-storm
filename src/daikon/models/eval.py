import numpy as np
import numpy.typing as npt
import torch
import torch.nn as nn
from sklearn.metrics import auc, confusion_matrix, f1_score, roc_curve
from sklearn.preprocessing import label_binarize
from torch.utils.data import DataLoader
from torchvision import transforms

from daikon.models import BrainTumorDataset


def test_model(
    model: nn.Module,
    transform: transforms.Compose,
    test_dir: str,
    device: torch.device,
    batch_size: int = 64,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], list[str]]:
    """
    Test the model.

    Args:
        model (nn.Module): Model to test. Model should already be on device
        transform (transforms.Compose): Data transformation to run
        test_dir (str): Directory containing test data
        device (torch.device): Device that the model is on
        batch_size (int): Number of images to test the model on at once

    Return:
        np.ndarray, np.ndarray: Predicted values, ground truth labels
    """
    # Get data
    test_dataset = BrainTumorDataset(test_dir, transform)
    test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    # Get output labels
    classes = test_dataset.get_classes()

    # Test the model
    pred = []
    ground_truth = []
    model.eval()
    with torch.no_grad():
        for inputs, labels in test_dataloader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)

            pred.append(outputs.detach().cpu().numpy())
            ground_truth.append(labels.detach().cpu().numpy())

    pred = np.vstack(pred)
    ground_truth = np.concatenate(ground_truth)
    return pred, ground_truth, classes


def get_confusion_f1(
    pred: npt.NDArray[np.float64],
    ground_truth: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], float]:
    """
    Get confusion matrix
    """
    pred_class = np.argmax(pred, axis=1)
    cm = confusion_matrix(ground_truth, pred_class)
    f1 = f1_score(ground_truth, pred_class, average="macro")

    return cm, f1


def get_roc(
    y_pred: npt.NDArray[np.float64],
    ground_truth: npt.NDArray[np.float64],
    classes: list[str],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float]:
    """
    Get the ROC curve and the area under the curve
    """
    y_true = label_binarize(ground_truth, classes=range(len(classes)))

    # Compute micro-average ROC
    fpr, tpr, _ = roc_curve(y_true.ravel(), y_pred.ravel())
    roc_auc = auc(fpr, tpr)

    return fpr, tpr, roc_auc
