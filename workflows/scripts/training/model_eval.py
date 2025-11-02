from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import torch

from daikon.models.dataclasses import ModelEvaluation, TrainingResults
from daikon.models.eval import get_confusion_f1, get_roc, test_model
from daikon.models.models import get_model


def main(
    model_path: str,
    output_path: str,
    training_res_path: str,
    test_dir: str,
    cfg: dict[str, Any],
):
    """
    Evaluate and visualize model performance.

    Args:
        model_path (str): Path to the saved model weights
        output_path (str): Output file path
        training_res_path (str): Training results path
        test_dir (str): Directory containing test data
        cfg (dict): Training config
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA version: {torch.version.cuda}")
    print("PyTorch CUDA available:", torch.cuda.is_available())
    print("GPU count:", torch.cuda.device_count())
    if torch.cuda.is_available():
        print("GPU name:", torch.cuda.get_device_name(0))

    # Load model weights and restore model
    model_state_dict = torch.load(model_path, weights_only=True)
    training_results = TrainingResults.load(training_res_path)

    model, transform = get_model(training_results.model_name)
    model.load_state_dict(model_state_dict)
    model = model.to(device)

    # Parameters
    batch_size = cfg["batch_size"]

    # Run eval
    y_pred, y_true, classes = test_model(model, transform, test_dir, device, batch_size)

    # Get confusion, f1 score, and roc
    confusion, f1 = get_confusion_f1(y_pred, y_true)
    fpr, tpr, auc = get_roc(y_pred, y_true, classes)

    # Save results
    eval_results = ModelEvaluation(
        model_name=training_results.model_name,
        y_pred=np.array(y_pred, dtype=np.float32),
        y_true=np.array(y_true, dtype=np.float32),
        confusion=np.array(confusion, dtype=np.float32),
        f1=f1,
        fpr=np.array(fpr, dtype=np.float32),
        tpr=np.array(tpr, dtype=np.float32),
        auc_roc=auc,
    )
    eval_results.save(output_path)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(
        model_path=snakemake.input[0],
        output_path=snakemake.output[0],
        training_res_path=snakemake.input[1],
        test_dir=snakemake.input[2],
        cfg=snakemake.config["training"],
    )
