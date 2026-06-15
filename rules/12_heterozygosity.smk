rule heterozygosity_filter:
    input:
        bcf          = os.path.join(RESULTS_DIR, "03_variants", "filtered.bcf"),
        males_list   = os.path.join(RESULTS_DIR, "03_variants", "males.txt"),
        females_list = os.path.join(RESULTS_DIR, "03_variants", "females.txt")
    output:
        done       = os.path.join(RESULTS_DIR, "12_hetero", "heterozygosity", "heterozygosity.done")
    params:
        sites      = os.path.join(RESULTS_DIR, "12_hetero", "heterozygosity", "passing_sites.txt"),
        vcf        = os.path.join(RESULTS_DIR, "12_hetero", "heterozygosity", "heterozygosity_filtered.vcf.gz"),
        male_het   = temp(os.path.join(RESULTS_DIR, "12_hetero", "heterozygosity", "male_het.vcf.gz")),
        female_hom = temp(os.path.join(RESULTS_DIR, "12_hetero", "heterozygosity", "female_hom.vcf.gz"))
    resources:
        cpus_per_task=1, mem_mb_per_cpu=8000, runtime=60
    log:
        os.path.join(RESULTS_DIR, "logs", "12_heterozygosity", "heterozygosity_filter.log")
    envmodules:
        "BCFtools/1.18-GCC-12.3.0"
    shell:
        r"""
        mkdir -p $(dirname {params.vcf}) $(dirname {log})
        N_MALES=$(wc -l < {input.males_list})
        N_FEMALES=$(wc -l < {input.females_list})
        AN_MALES=$(( N_MALES * 2 ))
        AN_FEMALES=$(( N_FEMALES * 2 ))
        echo "Males: $N_MALES | Females: $N_FEMALES" > {log}

        bcftools view -S {input.males_list} -v snps -m2 -M2 {input.bcf} -Ou 2>> {log} \
        | bcftools view -i "N_PASS(GT=\"het\")=${{N_MALES}} && F_MISSING=0 && AN=${{AN_MALES}}" \
            -Oz -o {params.male_het} 2>> {log}
        bcftools index {params.male_het} 2>> {log}

        bcftools view -S {input.females_list} -v snps -m2 -M2 {input.bcf} -Ou 2>> {log} \
        | bcftools view -i "N_PASS(GT=\"hom\")=${{N_FEMALES}} && F_MISSING=0 && AN=${{AN_FEMALES}}" \
            -Oz -o {params.female_hom} 2>> {log}
        bcftools index {params.female_hom} 2>> {log}

        ISEC_DIR=$(dirname {params.sites})/isec_tmp
        bcftools isec -n=2 -p $ISEC_DIR {params.male_het} {params.female_hom} 2>> {log}
        awk 'BEGIN{{OFS="\t"}}{{print $1,$2}}' $ISEC_DIR/sites.txt > {params.sites}

        TARGETS={params.sites}
        bgzip -c {params.sites} > $TARGETS
        tabix -s1 -b2 -e2 $TARGETS
        bcftools view -R $TARGETS -Oz -o {params.vcf} {input.bcf} 2>> {log}
        bcftools index {params.vcf} 2>> {log}

        rm -rf $ISEC_DIR $TARGETS ${{TARGETS}}.tbi

        touch {output.done}
        """