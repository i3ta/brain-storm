rule preprocess_images:
    input:
        "../data/raw/{folder}"
    output:
        directory("../data/processed/{folder}")
    script:
        "../scripts/preprocessing/preprocess.py"
