# =============================================================================
# 08 - GWAS ANALYSIS: GENOME-WIDE ASSOCIATION STUDY FOR SEX DETERMINATION
#
# Workflow:
# 1. Convert filtered BCF to PLINK ped/map (chromosomes only, SNPs, MAF filter)
# 2. Create PLINK binary bed with sex phenotype
# 3. Run GEMMA association test (linear model)
# =============================================================================


rule gwas_vcf:
    """
    Convert filtered BCF to PLINK ped/map format.
    - Keep only chromosomes matching the configured prefix
    - Remove indels, filter by MAF and missingness
    - Rename chromosomes to numeric IDs for PLINK compatibility
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf")
    output:
        ped = os.path.join(RESULTS_DIR, "08_GWAS", "gwas.ped"),
        map = os.path.join(RESULTS_DIR, "08_GWAS", "gwas.map"),
        chr_map = os.path.join(RESULTS_DIR, "08_GWAS", "chr_map.txt")
    params:
        base_name = os.path.join(RESULTS_DIR, "08_GWAS", "gwas"),
        chrom_prefix = CHROM_PREFIX
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "08_GWAS", "gwas_vcf.log")
    envmodules:
        "VCFtools/0.1.16-GCC-12.3.0",
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        """
        mkdir -p $(dirname {output.ped})
        mkdir -p $(dirname {log})

        # Build chromosome-to-integer map for PLINK compatibility
        bcftools view -H {input.vcf} 2> {log} | \
        cut -f1 | \
        grep '^{params.chrom_prefix}' | \
        sort -u | \
        awk '{{print $1 "\t" NR}}' > {output.chr_map} 2>> {log}

        # Rename chromosomes and convert to PLINK ped/map
        bcftools annotate --rename-chrs {output.chr_map} {input.vcf} 2>> {log} | \
        vcftools --vcf - \
            --plink \
            --remove-indels \
            --max-missing 0.5 \
            --max-maf 0.95 \
            --maf 0.05 \
            --out {params.base_name} 2>> {log}
        """


rule gwas_plink:
    """
    Create PLINK binary bed files with sex as phenotype.
    Males (matching male_pattern) = 1, females (matching female_pattern) = 2.
    """
    input:
        vcf = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf"),
        ped = os.path.join(RESULTS_DIR, "08_GWAS", "gwas.ped"),
        map = os.path.join(RESULTS_DIR, "08_GWAS", "gwas.map")
    output:
        bed = os.path.join(RESULTS_DIR, "08_GWAS", "gwas_plink.bed")
    params:
        base_name_in = os.path.join(RESULTS_DIR, "08_GWAS", "gwas"),
        base_name_out = os.path.join(RESULTS_DIR, "08_GWAS", "gwas_plink"),
        sex_meta = os.path.join(RESULTS_DIR, "08_GWAS", "sex_meta.txt"),
        male_pattern = config.get("male_pattern", "ML"),
        female_pattern = config.get("female_pattern", "FL")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "08_GWAS", "gwas_plink.log")
    envmodules:
        "PLINK/2.00a3.7-foss-2022a",
        "BCFtools/1.15.1-GCC-11.3.0"
    shell:
        """
        mkdir -p $(dirname {log})

        # Create phenotype file from sample names
        bcftools query -l {input.vcf} | \
        awk -v mp="{params.male_pattern}" -v fp="{params.female_pattern}" '
        BEGIN {{OFS="\t"}} {{
            sex = ($0 ~ mp) ? 1 : ($0 ~ fp) ? 2 : 0;
            print $0, $0, sex
        }}' > {params.sex_meta} 2> {log}

        plink --file {params.base_name_in} \
            --pheno {params.sex_meta} \
            --make-bed \
            --out {params.base_name_out} \
            --noweb \
            --allow-no-sex 2>> {log}
        """


rule gemma:
    """
    Run GEMMA linear model association test (lm 2) for sex phenotype.
    GEMMA writes output to a local ./output/ directory, so we cd first
    and move results to the expected paths.
    """
    input:
        bed = os.path.join(RESULTS_DIR, "08_GWAS", "gwas_plink.bed")
    output:
        assoc = os.path.join(RESULTS_DIR, "08_GWAS", "gwas_gemma.assoc.txt"),
        txt = os.path.join(RESULTS_DIR, "08_GWAS", "gwas_gemma.log.txt")
    params:
        base_name_in = os.path.join(RESULTS_DIR, "08_GWAS", "gwas_plink"),
        output_directory = os.path.join(RESULTS_DIR, "08_GWAS")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=40000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "08_GWAS", "gemma.log")
    shell:
        """
        mkdir -p $(dirname {log})

        cd {params.output_directory} && \
        /user/brussel/109/vsc10945/home/scratch/Snakemake/SexDetection/software/gemma-0.98.5-linux-static-AMD64 \
            -bfile {params.base_name_in} \
            -lm 4 \
            -o temp_output &> {log}

        mv {params.output_directory}/output/temp_output.assoc.txt {output.assoc}
        mv {params.output_directory}/output/temp_output.log.txt {output.txt}
        rm -rf {params.output_directory}/output
        """
