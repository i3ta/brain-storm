rule train_model:
    input:
        "../data/processed/{folder}/train/",
        "../data/processed/{folder}/cv/",
    output:
        "../data/processed/{folder}/models/{model_name}.pkl",
        "results/{folder}/training/{model_name}.pkl",
    script:
        "../scripts/training/training.py"

rule plot_training_loss:
    input:
        "results/{folder}/training/{model_name}.pkl"
    output:
        "results/{folder}/training/{model_name}_loss.png"
    script:
        "../scripts/training/training_loss.py"

rule run_model_eval:
    input:
        "../data/processed/{folder}/models/{model_name}.pkl",
        "results/{folder}/training/{model_name}.pkl",
        "../data/processed/{folder}/test/"
    output:
        "results/{folder}/eval/{model_name}_eval.pkl"
    script:
        "../scripts/training/model_eval.py"

rule run_model_eval_vis:
    input:
        "results/{folder}/eval/{model_name}_eval.pkl"
    output:
        "results/{folder}/eval/{model_name}_eval.png"
    script:
        "../scripts/training/model_eval_vis.py"
