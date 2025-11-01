from typing import Any

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torchvision.models as models
from torch.utils.data import DataLoader
from torchvision import transforms
from tqdm import tqdm

from daikon.models import BrainTumorDataset
from daikon.models.dataclasses import TrainingResults


def main(
    model_name: str,
    train_dir: str,
    cv_dir: str,
    model_path: str,
    training_data_path: str,
    cfg: dict[str, Any],
):
    """
    Train a model.

    Args:
        model (str): Name of the model to train
        train_dir (str): Training data
        cv_dir (str): Cross-validation data
        model_path (str): Path to save model to
        training_data (str): Path to save training data to
        cfg (dict): Training config
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA version: {torch.version.cuda}")
    print("PyTorch CUDA available:", torch.cuda.is_available())
    print("GPU count:", torch.cuda.device_count())
    if torch.cuda.is_available():
        print("GPU name:", torch.cuda.get_device_name(0))

    # TODO: Add each model here. The model just takes in a string, so you can pass
    #       parameters as part of the string to specify parameters for
    #       hyperparameter tuning.
    if model_name == "vgg16_pretrained":
        model = models.vgg16(weights=models.VGG16_Weights.IMAGENET1K_V1)
        model.classifier[6] = nn.Linear(model.classifier[6].in_features, 4)
        model = model.to(device)
        transform = transforms.Compose(
            [transforms.Resize((224, 224)), transforms.ToTensor()]
        )
    else:
        raise NotImplementedError(f"The model {model_name} cannot be found")

    # Training parameters
    epochs = cfg["epochs"]
    batch_size = cfg["batch_size"]
    learning_rate = cfg["learning_rate"]
    patience = cfg["patience"]

    # Get data
    train_dataset = BrainTumorDataset(train_dir, transform)
    train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    cv_dataset = BrainTumorDataset(cv_dir, transform)
    cv_dataloader = DataLoader(cv_dataset, batch_size=batch_size, shuffle=True)

    # Set up loss function and optimizer
    loss_fn = nn.CrossEntropyLoss().to(device)
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    # Run training loop
    train_loss = []
    best_val_loss = float("inf")
    epochs_no_impr = 0
    for _ in tqdm(range(epochs)):
        model.train()
        for data, labels in train_dataloader:
            data = data.to(device)
            labels = labels.to(device)

            # Predict with model and get loss
            labels_pred = model(data)
            loss = loss_fn(labels_pred, labels)

            # Backward pass and optimization
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        # Validate
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for data, labels in cv_dataloader:
                data = data.to(device)
                labels = labels.to(device)

                labels_pred = model(data)
                loss = loss_fn(labels_pred, labels)
                val_loss += loss.item()
        train_loss.append(val_loss)

        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            epochs_no_impr = 0
        else:
            epochs_no_impr += 1
            if epochs_no_impr >= patience:
                break

    # Save data
    training_results = TrainingResults(
        model_name=model_name,
        loss=np.array(train_loss, dtype=np.float32),
        batch_size=batch_size,
        learning_rate=learning_rate,
        patience=patience,
    )
    training_results.save(training_data_path)

    model = model.to("cpu")
    torch.save(model.state_dict(), model_path)


if __name__ == "__main__" and "snakemake" in locals():
    from snakemake.script import snakemake

    main(
        model_name=snakemake.wildcards.model_name,
        train_dir=snakemake.input[0],
        cv_dir=snakemake.input[1],
        model_path=snakemake.output[0],
        training_data_path=snakemake.output[1],
        cfg=snakemake.config["training"],
    )
