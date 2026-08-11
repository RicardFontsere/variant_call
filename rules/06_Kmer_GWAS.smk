# =============================================================================
# 06 - KMER GWAS ANALYSIS
#
# Workflow:
# 1. KMC3 k-mer counting (canonical + all) per sample
# 2. Add strand information per sample
# 3. Combine k-mers across samples
# 4. Build k-mer presence/absence table
# 5. Convert to PLINK bed format with phenotype
# 6. PLINK association test -> sex-specific k-mers
# 7. ABySS assembly of sex-specific k-mers
# 8. BLAST assembled contigs against reference
# =============================================================================

rule fastp_kmer:
    """
    Preprocess BGI/MGI reads for the k-mer branch.
    BGI/MGI libraries show non-random base composition in the first
    ~10-15 bp (random-hexamer priming bias), which produces spurious
    k-mers. fastp performs adapter + quality trimming and hard-crops the
    first 15 bp of both mates.
    """
    input:
        # ancient(): see rule fastp in 01_trimming.smk - raw FASTQs are too
        # large to checksum, so a bumped mtime alone would re-trigger the
        # whole k-mer branch.
        r1 = lambda wildcards: ancient(get_read_file(wildcards.sample, "1")),
        r2 = lambda wildcards: ancient(get_read_file(wildcards.sample, "2"))
    output:
        # temp(): deleted once rule kmer_analysis has consumed them. The QC
        # html/json below are small and kept as the permanent QC record.
        r1 = temp(os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1_kmer.fq.gz")),
        r2 = temp(os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2_kmer.fq.gz")),
        html = os.path.join(RESULTS_DIR, "00_qc", "{sample}_kmer.html"),
        json = os.path.join(RESULTS_DIR, "00_qc", "{sample}_kmer.json")
    params:
        crop_front = 15
    resources:
        cpus_per_task=10,
        mem_mb_per_cpu=2000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "fastp_kmer_{sample}.log")
    envmodules:
        "fastp/1.0.1-GCC-13.3.0"
    shell:
        """
        mkdir -p $(dirname {output.r1})
        mkdir -p $(dirname {output.json})
        mkdir -p $(dirname {log})
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --trim_front1 {params.crop_front} \
            --trim_front2 {params.crop_front} \
            --thread {resources.cpus_per_task} \
            -j {output.json} -h {output.html} &> {log}
        """

rule kmer_analysis:
    """
    Count k-mers (k=31) per sample using KMC3.
    Produces canonical (ci2) and all (ci0) k-mer databases.
    """
    input:
        r1 = os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_1_kmer.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "01_trimmed", "{sample}_2_kmer.fq.gz")
    output:
        canon = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_canon.kmc_pre"),
        all_kmers = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_all.kmc_pre")
    params:
        outdir = os.path.join(RESULTS_DIR, "06_kmer", "{sample}")
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=4000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "{sample}_kmer_analysis.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        cd {params.outdir}
        echo {input.r1} > input_files.txt
        echo {input.r2} >> input_files.txt
        {SOFTWARE_DIR}/kmc_v3 \
            -t{resources.cpus_per_task} -k31 -ci2 \
            @input_files.txt output_kmc_canon ./ \
            1> kmc_canon.1 2> kmc_canon.2
        {SOFTWARE_DIR}/kmc_v3 \
            -t{resources.cpus_per_task} -k31 -ci0 -b \
            @input_files.txt output_kmc_all ./ \
            1> kmc_all.1 2> kmc_all.2
        cat kmc_canon.2 kmc_all.2 > {log} 2>&1
        """


rule add_strand_information:
    """
    Add strand information to k-mers using canonical and all k-mer databases.
    """
    input:
        canon = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_canon.kmc_pre"),
        all_kmers = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_all.kmc_pre")
    output:
        strands = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "kmers_with_strand")
    params:
        canon_prefix = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_canon"),
        all_prefix = os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "output_kmc_all")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=100000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "{sample}_strand_info.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p $(dirname {log})
        {SOFTWARE_DIR}/kmers_add_strand_information \
            -c {params.canon_prefix} \
            -n {params.all_prefix} \
            -k 31 \
            -o {output.strands} &> {log}
        """


rule combine_kmers:
    """
    Combine k-mers found across multiple samples.
    Filters by minor allele count (mac=3) and presence proportion (p=0.2).
    """
    input:
        strands = expand(os.path.join(RESULTS_DIR, "06_kmer", "{sample}", "kmers_with_strand"), sample=SAMPLES)
    output:
        combined = os.path.join(RESULTS_DIR, "06_kmer", "combined", "intermediate", "combined_kmers"),
        kmers_list = os.path.join(RESULTS_DIR, "06_kmer", "combined", "intermediate", "kmers_to_combine")
    params:
        samples = SAMPLES
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "combine_kmers.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p $(dirname {output.kmers_list})
        mkdir -p $(dirname {log})
        paste <(printf '%s\\n' {input.strands}) <(printf '%s\\n' {params.samples}) > {output.kmers_list}
        {SOFTWARE_DIR}/list_kmers_found_in_multiple_samples \
            -l {output.kmers_list} \
            -k 31 \
            --mac 3 \
            -p 0.2 \
            -o {output.combined} &> {log}
        """


rule build_kmers_table:
    """
    Build k-mer presence/absence table across all samples.
    """
    input:
        combined = os.path.join(RESULTS_DIR, "06_kmer", "combined", "intermediate", "combined_kmers"),
        kmers_list = os.path.join(RESULTS_DIR, "06_kmer", "combined", "intermediate", "kmers_to_combine")
    output:
        table = os.path.join(RESULTS_DIR, "06_kmer", "combined", "tables", "kmers_table.table"),
        names = os.path.join(RESULTS_DIR, "06_kmer", "combined", "tables", "kmers_table.names")
    params:
        out_prefix = os.path.join(RESULTS_DIR, "06_kmer", "combined", "tables", "kmers_table")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000,
        runtime=1440
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "build_kmers_table.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p $(dirname {output.table})
        mkdir -p $(dirname {log})
        {SOFTWARE_DIR}/build_kmers_table \
            -l {input.kmers_list} \
            -k 31 \
            -a {input.combined} \
            -o {params.out_prefix} &> {log}
        """


rule kmers_table_to_bed:
    """
    Convert k-mer table to PLINK bed format.
    Creates phenotype file: males (matching male_pattern) = 1, females = 2.
    """
    input:
        table = os.path.join(RESULTS_DIR, "06_kmer", "combined", "tables", "kmers_table.table"),
        kmers_list = os.path.join(RESULTS_DIR, "06_kmer", "combined", "intermediate", "kmers_to_combine")
    output:
        bed = os.path.join(RESULTS_DIR, "06_kmer", "combined", "plink", "kmers_to_plink.0.bed")
    params:
        phenotype = os.path.join(RESULTS_DIR, "06_kmer", "combined", "intermediate", "phenotype.pheno"),
        table_prefix = os.path.join(RESULTS_DIR, "06_kmer", "combined", "tables", "kmers_table"),
        out_prefix = os.path.join(RESULTS_DIR, "06_kmer", "combined", "plink", "kmers_to_plink"),
        male_pattern = config.get("male_pattern", "ML")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=90000,
        runtime=2880
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "kmers_to_bed.log")
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p $(dirname {output.bed})
        mkdir -p $(dirname {log})
        echo -e "accession_id\\tphenotype_value" > {params.phenotype}
        awk -v mp="{params.male_pattern}" '{{print $2 "\\t" ($2~mp?1:2)}}' \
            {input.kmers_list} >> {params.phenotype}
        {SOFTWARE_DIR}/kmers_table_to_bed \
            -t {params.table_prefix} \
            -k 31 \
            -p {params.phenotype} \
            --maf 0.05 --mac 3 \
            -b 1000000000 \
            -o {params.out_prefix} &> {log}
        """


rule kmer_association:
    """
    Run PLINK association test on k-mer data.
    Extract strictly sex-specific k-mers:
        Male-specific:   present in ALL males (F_A=1), absent in ALL females (F_U=0)
        Female-specific: absent in ALL males (F_A=0), present in ALL females (F_U=1)

    PLINK .assoc has leading whitespace; we clean to TSV first so columns are:
        $1=CHR, $2=SNP, $3=BP, $4=A1, $5=F_A(males), $6=F_U(females),
        $7=A2, $8=CHISQ, $9=P, $10=OR
    """
    input:
        bed = os.path.join(RESULTS_DIR, "06_kmer", "combined", "plink", "kmers_to_plink.0.bed")
    output:
        assoc = os.path.join(RESULTS_DIR, "06_kmer", "combined", "plink", "kmers.assoc"),
        assoc_clean = os.path.join(RESULTS_DIR, "06_kmer", "combined", "plink", "kmers.assoc.tsv"),
        male_assoc = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "male_specific_kmers.assoc"),
        female_assoc = os.path.join(RESULTS_DIR, "06_kmer", "combined", "results", "female_specific_kmers.assoc"),
        male_abyss = os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "male_abyss.input"),
        female_abyss = os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "female_abyss.input")
    params:
        workdir = os.path.join(RESULTS_DIR, "06_kmer", "combined", "plink")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=32000,
        runtime=2000
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "kmer_association.log")
    envmodules:
        "PLINK/2.00a3.7-foss-2022a",
        "Python/3.10.4-GCCcore-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.male_assoc})
        mkdir -p $(dirname {output.male_abyss})
        mkdir -p $(dirname {log})
        cd {params.workdir}

        # Run association test
        plink --noweb --bfile kmers_to_plink.0 \
            --allow-no-sex --assoc --out kmers &> {log}

        # Clean leading whitespace to proper TSV
        sed -E 's/^[[:space:]]+//; s/[[:space:]]+/\\t/g' kmers.assoc > {output.assoc_clean} 2>> {log}

        # Male-specific: F_A=1 (all males have it), F_U=0 (no females have it)
        awk -F'\\t' 'NR>1 && $5 == 1 && $6 == 0' {output.assoc_clean} > {output.male_assoc} 2>> {log}
        # Female-specific: F_A=0 (no males have it), F_U=1 (all females have it)
        awk -F'\\t' 'NR>1 && $5 == 0 && $6 == 1' {output.assoc_clean} > {output.female_assoc} 2>> {log}

        # Log counts
        echo "Total associations: $(tail -n +2 {output.assoc_clean} | wc -l)" >> {log}
        echo "Male-specific (F_A=1, F_U=0): $(wc -l < {output.male_assoc})" >> {log}
        echo "Female-specific (F_A=0, F_U=1): $(wc -l < {output.female_assoc})" >> {log}

        # Convert to ABySS input
        python3 {SCRIPTS_DIR}/plink_to_abyss_kmers.py {output.male_assoc} {output.male_abyss} 2>> {log}
        python3 {SCRIPTS_DIR}/plink_to_abyss_kmers.py {output.female_assoc} {output.female_abyss} 2>> {log}
        """


rule abyss_male:
    """
    Assemble male-specific k-mers with ABySS.
    """
    input:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "male_abyss.input")
    output:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "male_abyss.output")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "abyss_male.log")
    envmodules:
        "ABySS/2.3.7-foss-2023a"
    shell:
        """
        mkdir -p $(dirname {log})
        if [ -s {input} ]; then
            ABYSS -k25 -c0 -e0 {input} -o {output} &> {log}
        else
            echo "Input {input} is empty - ...; writing empty output." > {log}
            touch {output}
        fi
        """


rule abyss_female:
    """
    Assemble female-specific k-mers with ABySS.
    """
    input:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "female_abyss.input")
    output:
        os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "female_abyss.output")
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "abyss_female.log")
    envmodules:
        "ABySS/2.3.7-foss-2023a"
    shell:
        """
        mkdir -p $(dirname {log})
        if [ -s {input} ]; then
            ABYSS -k25 -c0 -e0 {input} -o {output} &> {log}
        else
            echo "Input {input} is empty - ...; writing empty output." > {log}
            touch {output}
        fi
        """


rule blast_db:
    """
    Create BLAST nucleotide database from reference genome.
    """
    input:
        ref = REFERENCE
    output:
        db = REFERENCE + ".nin"
    resources:
        cpus_per_task=2,
        mem_mb_per_cpu=8000,
        runtime=120
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "blast_db.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        mkdir -p $(dirname {log})
        makeblastdb -in {input.ref} -dbtype nucl 2> {log}
        """


rule blast_male_kmers:
    """
    BLAST male-assembled contigs against reference genome.
    """
    input:
        contigs = os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "male_abyss.output"),
        ref = REFERENCE,
        db = REFERENCE + ".nin"
    output:
        blast = os.path.join(RESULTS_DIR, "06_kmer", "combined", "blast", "male_blast.out")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=15000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "blast_male.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        mkdir -p $(dirname {output.blast})
        mkdir -p $(dirname {log})
        if [ -s {input.contigs} ]; then
            blastn \
                -query {input.contigs} \
                -db {input.ref} \
                -task blastn \
                -perc_identity 80 \
                -qcov_hsp_perc 80 \
                -dust yes \
                -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
                -evalue 1e-10 \
                -num_threads {resources.cpus_per_task} \
                > {output.blast} 2> {log}
        else
            echo "Input {input.contigs} is empty - no contigs to BLAST; writing empty output." > {log}
            touch {output.blast}
        fi
        """


rule blast_female_kmers:
    """
    BLAST female-assembled contigs against reference genome.
    """
    input:
        contigs = os.path.join(RESULTS_DIR, "06_kmer", "combined", "assembly", "female_abyss.output"),
        ref = REFERENCE,
        db = REFERENCE + ".nin"
    output:
        blast = os.path.join(RESULTS_DIR, "06_kmer", "combined", "blast", "female_blast.out")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=4000,
        runtime=240
    log:
        os.path.join(RESULTS_DIR, "logs", "06_kmer", "blast_female.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        mkdir -p $(dirname {output.blast})
        mkdir -p $(dirname {log})
        if [ -s {input.contigs} ]; then
                    mkdir -p $(dirname {output.blast})
        mkdir -p $(dirname {log})
        if [ -s {input.contigs} ]; then
            blastn \
                -query {input.contigs} \
                -db {input.ref} \
                -task blastn \
                -perc_identity \
                -dust yes \
                -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
                -evalue 1e-10 \
                -num_threads {resources.cpus_per_task} \
                > {output.blast} 2> {log}
        else
            echo "Input {input.contigs} is empty - no contigs to BLAST; writing empty output." > {log}
            touch {output.blast}
        fi
        """