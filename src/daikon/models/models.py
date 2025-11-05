import torch.nn as nn
import torchvision.models as models
from torchvision import transforms


def get_model(model_name: str) -> tuple[nn.Module, transforms.Compose]:
    """
    Get the model and transform for the ML model.

    Args:
        model_name (str): Name of the model

    Return:
        model, transform
    """
    # TODO: Add each model here. The model just takes in a string, so you can pass
    #       parameters as part of the string to specify parameters for
    #       hyperparameter tuning.
    if model_name == "vgg16_pretrained":
        model = models.vgg16(weights=models.VGG16_Weights.IMAGENET1K_V1)
        model.classifier[6] = nn.Linear(model.classifier[6].in_features, 4)

        # Resize image to 224x224 for VGGNet
        transform = transforms.Compose(
            [
                transforms.Resize((224, 224)),
                transforms.ToTensor(),
                transforms.Normalize(
                    mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]
                ),
            ]
        )
    elif model_name == "vgg16":
        model = models.vgg16(weights=None)
        model.classifier[6] = nn.Linear(model.classifier[6].in_features, 4)
        transform = transforms.Compose(
            [
                transforms.Resize((224, 224)),
                transforms.ToTensor(),
                transforms.Normalize(
                    mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]
                ),
            ]
        )
    elif model_name == "resnet18_pretrained":
        model = models.resnet18(weights=models.ResNet18_Weights.IMAGENET1K_V1)
        model.fc = nn.Linear(model.fc.in_features, 4)
        
        transform = transforms.Compose(
            [
                transforms.ToTensor(),
                transforms.Normalize(
                    mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]
                ),
            ]
        )
    elif model_name == "resnet18":
        model = models.resnet18(weights=None)
        model.fc = nn.Linear(model.fc.in_features, 4)
        transform = transforms.Compose(
            [
                transforms.ToTensor()
            ]
        )
    else:
        raise NotImplementedError(f"The model {model_name} cannot be found")

    return model, transform
