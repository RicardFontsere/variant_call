"""
Variant calling and SNP analysis rules
"""

rule variant_calling:
    input:
        bams = expand(os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam"), sample=SAMPLES),
        bais = expand(os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam.csi"), sample=SAMPLES),
        ref = REFERENCE_GENOME
    output:
        vcf = expand(os.path.join(RESULTS_DIR, "variants", "variants_{species}.vcf.gz"), species=SPECIES)
    params:
        bamlist = os.path.join(RESULTS_DIR, "variants", "bamlist.txt"),
        min_base_qual = config["variant_calling"]["min_base_qual"],
        min_map_qual = config["variant_calling"]["min_map_qual"]
    resources:
        cpus_per_task = 64,
        mem_mb_per_cpu = 4000
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "variant_calling.log")
    envmodules:
        "BCFtools/1.21-GCC-13.3.0"
    shell:
        """
        mkdir -p {RESULTS_DIR}/variants
        mkdir -p {RESULTS_DIR}/logs/variants
        ls -1 {RESULTS_DIR}/mapped/*.sorted.bam > {params.bamlist}
        bcftools mpileup --threads {resources.cpus_per_task} -f {input.ref} \
          --min-BQ {params.min_base_qual} \
          --min-MQ {params.min_map_qual} \
          --annotate FORMAT/DP,FORMAT/AD,FORMAT/SP \
          --bam-list {params.bamlist} | \
        bcftools call --threads {resources.cpus_per_task} -mv -O z -o {output.vcf} 2> {log}
        """

rule index_vcf:
    input:
        vcf = os.path.join(RESULTS_DIR, "variants", "variants_{prefix}.vcf.gz")
    output:
        index = os.path.join(RESULTS_DIR, "variants", "variants_{prefix}.vcf.gz.tbi")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "index_{prefix}.log")
    envmodules:
        "BCFtools/1.21-GCC-13.3.0"
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=1000
    shell:
        """
        bcftools index -c {input.vcf} -o {output.index} 2> {log}
        """

rule filter_variants:
    input:
        vcf = os.path.join(RESULTS_DIR, "variants", "variants_{species}.vcf.gz"),
        index = os.path.join(RESULTS_DIR, "variants", "variants_{species}.vcf.gz.tbi")
    output:
        vcf = os.path.join(RESULTS_DIR, "variants", "variants_filtered_{species}.vcf.gz")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "variant_filtering_{species}.log")
    envmodules:
        "BCFtools/1.21-GCC-13.3.0"
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000
    shell:
        """
        bcftools view {input.vcf} --include 'QUAL >= 30 && TYPE="snp" && N_ALT=1 && INFO/DP>=60 && AC>=2 && COUNT(GT="alt")>=2 && COUNT(FMT/DP>=10)=6 && COUNT(FMT/DP<=40)=6' -o {output.vcf}
        """

rule split_vcfs:
    input:
        filtered_vcf = expand(os.path.join(RESULTS_DIR, "variants", "variants_filtered_{species}.vcf.gz"), species=SPECIES)
    output:
        sample_vcfs = expand(os.path.join(RESULTS_DIR, "variants", "samples", "{sample}.vcf.gz"), sample=SAMPLES)
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "variant_splitting.log")
    params:
        samples_dir = os.path.join(RESULTS_DIR, "variants", "samples")
    envmodules:
        "BCFtools/1.21-GCC-13.3.0"
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000
    shell:
        """
        mkdir -p {params.samples_dir}
        samples=$(bcftools query -l {input.filtered_vcf})
        
        for full_path in $samples; do
            bam_file=$(basename "$full_path")
            sample=$(echo "$bam_file" | sed 's/\.sorted\.bam$//')
            echo "Extracting sample $sample from VCF (from $full_path)" >> {log}
            bcftools view -c1 -s "$full_path" \
            --types snps \
            -O z \
            -o "{params.samples_dir}/${{sample}}.vcf.gz" \
            {input.filtered_vcf}
        done
        """

rule SNP_density:
    input:
        split_vcf = os.path.join(RESULTS_DIR, "variants", "samples", "{sample}.vcf.gz"),
        windows = WINDOWS_FILE
    output:
        snp_density = os.path.join(RESULTS_DIR, "variants", "samples", "{sample}_50kb.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "snp_density_{sample}.log")
    envmodules:
        "BEDTools/2.30.0-GCC-11.3.0"
    resources:
        cpus_per_task=10,
        mem_mb_per_cpu=4000
    shell:
        """
        mkdir -p $(dirname {output.snp_density})
        mkdir -p $(dirname {log})
        bedtools coverage \
            -a {input.windows} \
            -b {input.split_vcf} > {output.snp_density} 2> {log}
        """

rule normalize_snp_male:
    input:
        snp_density = os.path.join(RESULTS_DIR, "variants", "samples", "{sample}_50kb.txt")
    output:
        norm = os.path.join(RESULTS_DIR, "variants", "samples", "malenormalized", "{sample}_50kb_normalized.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "snp_normalize_male_{sample}.log")
    envmodules:
        "R/4.2.1-foss-2022a"
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000
    shell:
        """
        mkdir -p $(dirname {output.norm})
        Rscript --vanilla -e '
            data <- read.table("{input.snp_density}", sep="\\t", header=FALSE);
            data_sorted <- data[order(data$V1, data$V2), ];
            data_sorted$V5 <- data_sorted$V4 / median(data_sorted$V4);
            write.table(data_sorted, "{output.norm}", sep="\\t", 
                      row.names=FALSE, col.names=FALSE, quote=FALSE)
        ' &> {log}
        """

rule normalize_snp_female:
    input:
        snp_density = os.path.join(RESULTS_DIR, "variants", "samples", "{sample}_50kb.txt")
    output:
        norm = os.path.join(RESULTS_DIR, "variants", "samples", "femalenormalized", "{sample}_50kb_normalized.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "snp_normalize_female_{sample}.log")
    envmodules:
        "R/4.2.1-foss-2022a"
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000
    shell:
        """
        mkdir -p $(dirname {output.norm})
        Rscript --vanilla -e '
            data <- read.table("{input.snp_density}", sep="\\t", header=FALSE);
            data_sorted <- data[order(data$V1, data$V2), ];
            data_sorted$V5 <- data_sorted$V4 / median(data_sorted$V4);
            write.table(data_sorted, "{output.norm}", sep="\\t", 
                      row.names=FALSE, col.names=FALSE, quote=FALSE)
        ' &> {log}
        """

rule calculate_snp_density_male:
    input:
        normalized_files = expand(os.path.join(RESULTS_DIR, "variants", "samples", "malenormalized", "{sample}_50kb_normalized.txt"), sample=ML_SAMPLES)
    output:
        male_snp = os.path.join(RESULTS_DIR, "variants", "snpdensity_males.txt")
    params:
        male_dir = os.path.join(RESULTS_DIR, "variants", "samples", "malenormalized")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "snp_density_males.log")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000
    envmodules:
        "Python/3.10.4-GCCcore-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.male_snp})
        python3 /user/brussel/109/vsc10945/home/scratch/Snakemake/SexDetection/scripts/extract_snp_density.py \
            {params.male_dir} {output.male_snp} \
            &> {log}
        """

rule calculate_snp_density_female:
    input:
        normalized_files = expand(os.path.join(RESULTS_DIR, "variants", "samples", "femalenormalized", "{sample}_50kb_normalized.txt"), sample=FL_SAMPLES)
    output:
        female_snp = os.path.join(RESULTS_DIR, "variants", "snpdensity_females.txt")
    params:
        female_dir = os.path.join(RESULTS_DIR, "variants", "samples", "femalenormalized")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "snp_density_females.log")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000
    envmodules:
        "Python/3.10.4-GCCcore-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.female_snp})
        python3 /user/brussel/109/vsc10945/home/scratch/Snakemake/SexDetection/scripts/extract_snp_density.py \
            {params.female_dir} {output.female_snp} \
            &> {log}
        """

rule calculate_snp_density_fold_change:
    input:
        female_snp = os.path.join(RESULTS_DIR, "variants", "snpdensity_females.txt"),
        male_snp = os.path.join(RESULTS_DIR, "variants", "snpdensity_males.txt")
    output:
        snp_fc = os.path.join(RESULTS_DIR, "variants", "snpdensity_fc.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "variants", "snp_density_fc.log")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000
    envmodules:
        "Python/3.10.4-GCCcore-11.3.0"
    shell:
        """
        python3 /user/brussel/109/vsc10945/home/scratch/Snakemake/SexDetection/scripts/snpdensity_fold_change.py \
        {input.female_snp} \
        {input.male_snp} \
        {output.snp_fc} &> {log}
        """
