rule run_eda_folder:
    input:
        "../data/{folder}/{type}"
    output:
        "results/{folder}/summary/{type}_data_summary.pkl"
    script:
        "../scripts/eda_folder.py"

rule run_eda:
    input:
        "results/{folder}/summary/glioma_data_summary.pkl",
        "results/{folder}/summary/meningioma_data_summary.pkl",
        "results/{folder}/summary/no_tumor_data_summary.pkl",
        "results/{folder}/summary/pituitary_data_summary.pkl",
    output:
        "results/{folder}/summary/raw_data_summary.png"
    script:
        "../scripts/eda.py"

rule image_dim_dist:
    input:
        "results/{folder}/summary/{type}_data_summary.pkl",
    output:
        "results/{folder}/summary/{type}_dims.png"
    script:
        "../scripts/image_dim_dist.py"
