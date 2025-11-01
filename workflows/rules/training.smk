rule train_model:
    input:
        "../data/processed/{folder}/train/",
        "../data/processed/{folder}/cv/",
    output:
        "../data/{folder}/models/{model_name}.pkl",
        "results/{folder}/training/{model_name}.pkl",
    script:
        "../scripts/training/training.py"
