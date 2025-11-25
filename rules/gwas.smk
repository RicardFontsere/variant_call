"""
GWAS analysis rules: genome-wide association study for sex determination
"""

rule GWAS_vcf:
    input:
        vcf = os.path.join(RESULTS_DIR, "variants", "variants_{species}.vcf.gz")
    output:
        ped = os.path.join(RESULTS_DIR, "GWAS", "GWAS_{species}.ped"),
        map = os.path.join(RESULTS_DIR, "GWAS", "GWAS_{species}.map"),
        chr_map = os.path.join(RESULTS_DIR, "GWAS", "chr_map_{species}.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "gwas", "GWAS_filter_{species}.log")
    resources:
        cpus_per_task=1,
        mem_mb_per_cpu=4000
    envmodules:
        "VCFtools/0.1.16-GCC-12.3.0",
        "BCFtools/1.18-GCC-12.3.0"
    params:
        base_name = os.path.join(RESULTS_DIR, "GWAS", "GWAS_{species}"),
        prefix = config["INPUT"]["chr_prefix"]  # Prefix for chromosomes to include (e.g., "OZ" for OZ*)
    shell:
        """
        mkdir -p $(dirname {output.ped})
        mkdir -p $(dirname {log})
        
        # Extract and filter chromosomes starting with OZ
        bcftools view -H {input.vcf} 2>> {log} | \
        cut -f1 | \
        grep '^{params.prefix}' | \
        sort -u | \
        awk '{{print $1 "\t" NR}}' > {output.chr_map} 2>> {log}
        
        # Rest of the pipeline remains unchanged
        bcftools annotate --rename-chrs {output.chr_map} {input.vcf} 2>> {log} | \
        vcftools --vcf - \
            --plink \
            --remove-indels \
            --max-missing 0.5 \
            --max-maf 0.95 \
            --maf 0.05 \
            --out {params.base_name} 2>> {log}
        """

rule GWAS_plink:
    input:
        vcf = os.path.join(RESULTS_DIR, "variants", "variants_{species}.vcf.gz"),
        ped = os.path.join(RESULTS_DIR, "GWAS", "GWAS_{species}.ped"),
        map = os.path.join(RESULTS_DIR, "GWAS", "GWAS_{species}.map")   
    output:
        os.path.join(RESULTS_DIR, "GWAS", "GWAS_plink_{species}.bed")
    params:
        base_name_in = os.path.join(RESULTS_DIR, "GWAS", "GWAS_{species}"),
        base_name_out = os.path.join(RESULTS_DIR, "GWAS", "GWAS_plink_{species}"),
        sex_meta = os.path.join(RESULTS_DIR, "GWAS", "sex_meta_{species}.txt")
    log: 
        os.path.join(RESULTS_DIR, "logs", "gwas", "GWAS_plink_{species}.log")
    envmodules:
        "PLINK/2.00a3.7-foss-2022a",
        "BCFtools/1.15.1-GCC-11.3.0"
    shell:
        """
        bcftools query -l {input.vcf} | \
        awk 'BEGIN {{OFS="\t"}} {{
            # Check for ML (male=1) or FL (female=2)
            sex = ($0 ~ /ML/) ? 1 : ($0 ~ /FL/) ? 2 : 0;
            print $0, $0, sex
        }}' > {params.sex_meta}
        plink --file {params.base_name_in} --pheno {params.sex_meta} --make-bed --out {params.base_name_out} --noweb --allow-no-sex
        """

rule GEMMA:
    input:
        os.path.join(RESULTS_DIR, "GWAS", "GWAS_plink_{species}.bed")
    output:
        assoc = os.path.join(RESULTS_DIR, "GWAS", "GWAS_gemma_{species}.assoc.txt"),
        txt = os.path.join(RESULTS_DIR, "GWAS", "GWAS_gemma_{species}.log.txt")
    params:
        base_name_in = os.path.join(RESULTS_DIR, "GWAS", "GWAS_plink_{species}"),
        base_name_out = os.path.join(RESULTS_DIR, "GWAS", "GWAS_gemma_{species}"),
        output_directory = os.path.join(RESULTS_DIR, "GWAS")
    log:
        os.path.join(RESULTS_DIR, "logs", "gwas", "GWAS_gemma_{species}.log")
    resources:
        partition="skylake_mpi",
        cpus_per_task=4,
        mem_mb_per_cpu=40000
    envmodules:
        "legacy-software/until-2021b",
        "GCC/5.4.0-2.26"
    shell:
        """
        cd {params.output_directory} && 
        /user/brussel/109/vsc10945/home/scratch/Software/kmerGWAS/external_programs/gemma_0_96 -bfile {params.base_name_in} -lm 2 -o temp_output &> {log} 
        mv ./output/temp_output.assoc.txt {output.assoc}
        mv ./output/temp_output.log.txt {output.txt}
        rm -r output
        """