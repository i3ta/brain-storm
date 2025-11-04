rule preprocess_images:
    input:
        "../data/raw/{folder}"
    output:
        directory("../data/processed/{folder}/test/"),
        directory("../data/processed/{folder}/cv/"),
        directory("../data/processed/{folder}/train/"),
    script:
        "../scripts/preprocessing/preprocess.py"
