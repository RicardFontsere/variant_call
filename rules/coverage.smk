rule generate_genome_length:
    input:
        ref = REFERENCE_GENOME
    output:
        length = GENOME_LENGTH 
    log:
        os.path.join(RESULTS_DIR, "logs", "preparatory", "genome_length.log")
    shell:
        """
        awk '
            BEGIN {{ OFS="\t" }}
            /^>/ {{
                if (id && len > 0) print id, len
                id = substr($1, 2)
                len = 0
                next
            }}
            {{ len += length }}
            END {{ if (id && len > 0) print id, len }}
        ' {input.ref} > {output.length}
        echo "Generated genome length file: {output.length}" >> {log}
        """

rule create_windows:
    input:
        genome_length = GENOME_LENGTH
    output:
        windows = WINDOWS_FILE
    params:
        window = config["coverage_params"]["window_size"]
    log:
        os.path.join(RESULTS_DIR, "logs", "coverage", "create_windows.log")
    envmodules:
        "BEDTools/2.30.0-GCC-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.windows})
        bedtools makewindows -g {input.genome_length} -w {params.window} > {output.windows}
        echo "Generated genome length file: {output.windows}" >> {log}
        """

rule calculate_coverage:
    input:
        bam = os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam"),
        bai = os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam.csi"),
        windows = WINDOWS_FILE
    output:
        cov = os.path.join(RESULTS_DIR, "coverage", "{sample}_multicov.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "coverage", "{sample}.log")
    envmodules:
        "BEDTools/2.30.0-GCC-11.3.0"
    shell:
        """
        mkdir -p "$(dirname {output.cov})"
        mkdir -p "$(dirname {log})"
    
        bedtools multicov \
        -bams {input.bam} \
        -bed {input.windows} \
        > {output.cov} 2> {log}
        """

rule normalize_coverage_male:
    input:
        cov = os.path.join(RESULTS_DIR, "coverage", "{sample}_multicov.txt")
    output:
        norm = os.path.join(RESULTS_DIR, "coverage", "malenormalized", "{sample}_multicov_normalized.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "coverage", "{sample}_normalize.log")
    envmodules:
        "R/4.2.1-foss-2022a"
    shell:
        """
        mkdir -p $(dirname {output.norm})
        Rscript --vanilla -e '
            data <- read.table("{input.cov}", sep="\\t", header=FALSE);
            data_sorted <- data[order(data$V1, data$V2), ];
            data_sorted$V5 <- data_sorted$V4 / median(data_sorted$V4);
            write.table(data_sorted, "{output.norm}", sep="\\t", 
                      row.names=FALSE, col.names=FALSE, quote=FALSE)
        ' &> {log}
        """

rule normalize_coverage_female:
    input:
        cov = os.path.join(RESULTS_DIR, "coverage", "{sample}_multicov.txt")
    output:
        norm = os.path.join(RESULTS_DIR, "coverage", "femalenormalized", "{sample}_multicov_normalized.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "coverage", "{sample}_normalize.log")
    envmodules:
        "R/4.2.1-foss-2022a"
    shell:
        """
        mkdir -p $(dirname {output.norm})
        Rscript --vanilla -e '
            data <- read.table("{input.cov}", sep="\\t", header=FALSE);
            data_sorted <- data[order(data$V1, data$V2), ];
            data_sorted$V5 <- data_sorted$V4 / median(data_sorted$V4);
            write.table(data_sorted, "{output.norm}", sep="\\t", 
                      row.names=FALSE, col.names=FALSE, quote=FALSE)
        ' &> {log}
        """

rule calculate_sex_coverage_male:
    input:
        normalized_files = expand(os.path.join(RESULTS_DIR, "coverage", "malenormalized", "{sample}_multicov_normalized.txt"), sample=ML_SAMPLES)
    output:
        male_cov = os.path.join(RESULTS_DIR, "coverage", "coverage_males.txt")
    params:
        male_dir = os.path.join(RESULTS_DIR, "coverage", "malenormalized")
    log:
        os.path.join(RESULTS_DIR, "logs", "coverage_males.log")
    resources:
        cpus_per_task=1,
        time="00:10:00",
        mem_mb_per_cpu=4000
    envmodules:
        "Python/3.10.4-GCCcore-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.male_cov})
        python3 /user/brussel/109/vsc10945/home/scratch/Snakemake/SexDetection/scripts/extract_coverage.py \
            {params.male_dir} {output.male_cov} \
            &> {log}
        """

rule calculate_sex_coverage_female:
    input:
        normalized_files = expand(os.path.join(RESULTS_DIR, "coverage", "femalenormalized", "{sample}_multicov_normalized.txt"), sample=FL_SAMPLES)
    output:
        female_cov = os.path.join(RESULTS_DIR, "coverage", "coverage_females.txt")
    params:
        female_dir = os.path.join(RESULTS_DIR, "coverage", "femalenormalized")
    log:
        os.path.join(RESULTS_DIR, "logs", "coverage_females.log")
    envmodules:
        "Python/3.10.4-GCCcore-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.female_cov})
        python3 /user/brussel/109/vsc10945/home/scratch/Snakemake/SexDetection/scripts/extract_coverage.py \
            {params.female_dir} {output.female_cov} \
            &> {log}
        """

rule calculate_coverage_fold_change:
    input:
        female_cov = os.path.join(RESULTS_DIR, "coverage", "coverage_females.txt"),
        male_cov = os.path.join(RESULTS_DIR, "coverage", "coverage_males.txt")
    output:
        m2f_cov = os.path.join(RESULTS_DIR, "coverage","male2femalecov.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "m2f.log")
    envmodules:
        "Python/3.10.4-GCCcore-11.3.0"
    shell:
        """
        python3 /user/brussel/109/vsc10945/home/scratch/Snakemake/SexDetection/scripts/coverage_fold_change.py \
        {input.female_cov} \
        {input.male_cov} \
        {output.m2f_cov} &> {log}
        """