import matplotlib.pyplot as plt
import numpy as np

from daikon.models.dataclasses import TrainingResults


def main(results_path: str, output: str):
    """
    Graph the training loss for a specific model for a specific model.

    Args:
        results_path (str): Path to the training results pickle file
        output (str): Path to the output file
    """
    results = TrainingResults.load(results_path)

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    x_labels = np.arange(0, results.loss.shape[0])
    ax.plot(x_labels, results.loss)
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Loss")
    ax.set_title(
        f"Training loss vs. Epoch\n(Model: {results.model_name}; Epochs: {results.loss.shape[0]})",
    )

    fig.tight_layout()
    fig.savefig(output)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(results_path=snakemake.input[0], output=snakemake.output[0])
