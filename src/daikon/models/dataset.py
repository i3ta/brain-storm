import os

import torch
from PIL import Image
from torch.utils.data import Dataset

image_extensions = (".png", ".jpg", ".jpeg")


class BrainTumorDataset(Dataset):
    """
    Generalized class for dataset.
    """

    def __init__(self, data_dir: str, transform=None):
        super().__init__()
        self.data_dir = data_dir
        self.transform = transform

        # Get folders/classes - only include directories
        self.classes = [
            c
            for c in os.listdir(self.data_dir)
            if os.path.isdir(os.path.join(self.data_dir, c))
        ]

        # Get image paths and labels
        self.images = []
        self.labels = []
        for c in self.classes:
            src_dir = os.path.join(self.data_dir, c)
            if not os.path.isdir(src_dir):
                continue
            for image in os.listdir(src_dir):
                if image.lower().endswith(image_extensions):
                    # is an image
                    self.images.append(os.path.join(src_dir, image))
                    self.labels.append(self.get_class_index(c))

    def __len__(self):
        return len(self.images)

    def __getitem__(self, idx: int):
        # Load image
        image_path = self.images[idx]
        image = Image.open(image_path).convert("RGB")

        if self.transform:
            image = self.transform(image)

        label = torch.tensor(self.labels[idx], dtype=torch.long)
        return image, label

    def get_class_index(self, label: str) -> int:
        return self.classes.index(label)

    def get_class_name(self, index: int) -> str:
        return self.classes[index]

    def get_classes(self) -> list[str]:
        return self.classes
