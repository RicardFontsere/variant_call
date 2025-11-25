rule trimmomatic:
    input:
        r1 = lambda wildcards: get_resource_path("reads", wildcards.sample, "1"),
        r2 = lambda wildcards: get_resource_path("reads", wildcards.sample, "2")
    output:
        r1 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_1_trimmed.fq.gz"),
        r1_unpaired = os.path.join(RESULTS_DIR, "trimmed", "unpaired", "{sample}_1_unpaired.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_2_trimmed.fq.gz"),
        r2_unpaired = os.path.join(RESULTS_DIR, "trimmed", "unpaired", "{sample}_2_unpaired.fq.gz")
    resources:
        cpus_per_task = 64
    params:
        adapters = config["trimming_params"]["adapters"],
        slidingwindow = config["trimming_params"]["slidingwindow"],
        minlen = config["trimming_params"]["minlen"],
        leading = config["trimming_params"]["leading"],
        trailing = config["trimming_params"]["trailing"]
    log:
        os.path.join(RESULTS_DIR, "logs", "trimmomatic", "{sample}.log")
    envmodules:
        "Trimmomatic/0.39-Java-11"
    shell:
        """
        mkdir -p {RESULTS_DIR}/trimmed
        mkdir -p {RESULTS_DIR}/trimmed/unpaired
        mkdir -p {RESULTS_DIR}/logs/trimmomatic
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads {resources.cpus_per_task} -phred33 \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            SLIDINGWINDOW:{params.slidingwindow} \
            LEADING:{params.leading} \
            TRAILING:{params.trailing} \
            MINLEN:{params.minlen} &> {log}
        """
