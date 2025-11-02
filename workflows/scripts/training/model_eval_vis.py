import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import ConfusionMatrixDisplay

from daikon.models.dataclasses import ModelEvaluation


def main(model_eval_path: str, output_path: str):
    """
    Visualize the metrics from model evaluation.

    Args:
        model_eval_path (str): Model evaluation pickle file location
        output_path (str): Output figure path
    """
    # Load data
    model_eval: ModelEvaluation = ModelEvaluation.load(model_eval_path)

    # Set up figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

    # Plot confusion matrix
    test_n = np.sum(model_eval.confusion, axis=1, dtype=np.int32)[0]
    cmd = ConfusionMatrixDisplay(
        confusion_matrix=model_eval.confusion,
        display_labels=model_eval.labels,
    )
    cmd.plot(ax=ax1, cmap="Blues")
    ax1.set_title(f"Confusion Matrix (n = {test_n})")

    # Plot ROC curve
    ax2.plot(
        model_eval.fpr, model_eval.tpr, "b", label=f"AUC = {model_eval.auc_roc:.5f}"
    )
    ax2.legend()
    ax2.set_title("Receiver operating characteristic curve")
    ax2.set_xlabel("False Positive Rate")
    ax2.set_ylabel("True Positive Rate")

    # Save figure
    fig.suptitle(
        f"Model Evaluation Metrics\n(Model: {model_eval.model_name}; F-score: {model_eval.f1:.5f})"
    )
    fig.tight_layout()
    fig.savefig(output_path)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(model_eval_path=snakemake.input[0], output_path=snakemake.output[0])
