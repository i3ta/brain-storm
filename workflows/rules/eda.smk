CLASSES = ["glioma", "meningioma", "no_tumor", "pituitary"]
CLUSTERS = list(range(1, 19))

rule run_eda_folder:
    input:
        "../data/{folder}/{type}"
    output:
        "results/{folder}/summary/{type}_data_summary.pkl"
    script:
        "../scripts/eda/eda_folder.py"

rule run_eda:
    input:
        lambda wildcards: expand("results/{folder}/summary/{classes}_data_summary.pkl", folder=wildcards.folder, classes=CLASSES)
    output:
        "results/{folder}/summary/raw_data_summary.png"
    script:
        "../scripts/eda/eda.py"

rule run_image_dim_dist:
    input:
        lambda wildcards: expand("results/{folder}/summary/{classes}_data_summary.pkl", folder=wildcards.folder, classes=CLASSES)
    output:
        "results/{folder}/summary/raw_image_dims.png"
    script:
        "../scripts/eda/image_dim_dist.py"

rule run_class_histograms:
    input:
        "../data/{folder}/{type}"
    output:
        "results/{folder}/hist/data/{type}_hist.pkl"
    script:
        "../scripts/eda/class_histogram.py"

rule run_histogram_clustering:
    input:
        lambda wildcards: expand("results/{folder}/hist/data/{classes}_hist.pkl", folder=wildcards.folder, classes=CLASSES)
    output:
        "results/{folder}/hist/cluster/cluster_{k}.pkl"
    script:
        "../scripts/eda/histogram_clustering.py"

rule run_class_rings:
    input:
        "../data/{folder}/{type}"
    output:
        "results/{folder}/rings/data/{type}_hist.pkl"
    script:
        "../scripts/eda/class_rings.py"

rule run_rings_clustering:
    input:
        lambda wildcards: expand("results/{folder}/rings/data/{classes}_hist.pkl", folder=wildcards.folder, classes=CLASSES)
    output:
        "results/{folder}/rings/cluster/cluster_{k}.pkl"
    script:
        "../scripts/eda/rings_clustering.py"

rule run_visualize_clustering:
    input:
        lambda wildcards: expand("results/{folder}/{cluster_type}/cluster/cluster_{clusters}.pkl", 
          folder=wildcards.folder, cluster_type=wildcards.cluster_type, clusters=CLUSTERS),
        lambda wildcards: expand("../data/{folder}/{classes}", folder=wildcards.folder, classes=CLASSES)
    output:
        "results/{folder}/{cluster_type}/cluster_{k}.png",
        "results/{folder}/{cluster_type}/example_{k}.png"
    script:
        "../scripts/eda/visualize_clustering.py"
