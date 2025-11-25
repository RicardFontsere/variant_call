rule kmer_analysis:
    input:
        r1 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_1_trimmed.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_2_trimmed.fq.gz")
    output:
        os.path.join(RESULTS_DIR, "kmer", "{sample}", "output_kmc_canon.kmc_pre"),
        os.path.join(RESULTS_DIR, "kmer", "{sample}", "output_kmc_all.kmc_pre")
    params:
        outdir = os.path.join(RESULTS_DIR, "kmer", "{sample}")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "{sample}_kmer_analysis.log")
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=4000
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})
        cd {params.outdir}
        echo {input.r1} > input_files.txt
        echo {input.r2} >> input_files.txt
        /user/brussel/109/vsc10945/home/scratch/Software/kmerGWAS/external_programs/kmc_v3 -t20 -k31 -ci2 @input_files.txt output_kmc_canon ./ 1> kmc_canon.1 2> kmc_canon.2
        /user/brussel/109/vsc10945/home/scratch/Software/kmerGWAS/external_programs/kmc_v3 -t20 -k31 -ci0 -b @input_files.txt output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2 
        """

rule add_strand_information:
    input:
        canon = os.path.join(RESULTS_DIR, "kmer", "{sample}", "output_kmc_canon.kmc_pre"),
        all_kmers = os.path.join(RESULTS_DIR, "kmer", "{sample}", "output_kmc_all.kmc_pre")
    output:
        os.path.join(RESULTS_DIR, "kmer", "{sample}", "kmers_with_strand")
    params:
        outdir = os.path.join(RESULTS_DIR, "kmer", "{sample}"),
        canon_prefix = os.path.join(RESULTS_DIR, "kmer", "{sample}", "output_kmc_canon"),
        all_prefix = os.path.join(RESULTS_DIR, "kmer", "{sample}", "output_kmc_all"),
        outfile = os.path.join(RESULTS_DIR, "kmer", "{sample}", "kmers_with_strand")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "{sample}_kmerstrandinfo.log")
    resources:
        cpus_per_task=4,
        partition="zen4,zen5_mpi",
        mem_mb_per_cpu=100000
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        cd {params.outdir}
        /user/brussel/109/vsc10945/home/scratch/Software/kmerGWAS/bin/kmers_add_strand_information -c {params.canon_prefix} -n {params.all_prefix} -k 31 -o {params.outfile}  &> {log}
        """

rule combine_kmersGWAS:
    input:
        strands = expand(os.path.join(RESULTS_DIR, "kmer", "{sample}", "kmers_with_strand"), sample=SAMPLES)
    output:
        combined_kmers = os.path.join(RESULTS_DIR, "kmer", "combined", "combined_kmers"),
        kmers_list = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_to_combine")
    params:
        samples = SAMPLES 
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "combine_kmerGWAS.log")
    resources:
        cpus_per_task=4,
        partition="zen4,zen5_mpi",
        mem_mb_per_cpu=4000
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        mkdir -p $(dirname {output.kmers_list})
        mkdir -p $(dirname {log})
        paste <(printf '%s\\n' {input.strands}) <(printf '%s\\n' {params.samples}) > {output.kmers_list}
        /user/brussel/109/vsc10945/home/scratch/Software/kmerGWAS/bin/list_kmers_found_in_multiple_samples \
            -l {output.kmers_list} \
            -k 31 \
            --mac 3 \
            -p 0.2 \
            -o {output.combined_kmers} &> {log}
        """

rule table_kmersGWAS:
    input:
        combined_kmers = os.path.join(RESULTS_DIR, "kmer", "combined", "combined_kmers"),
        kmers_list = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_to_combine")
    output:
        os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table.table"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table.names")
    params:
        outfile = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "kmers_table.log")
    resources:
        cpus_per_task=4,
        partition="zen4,zen5_mpi",
        mem_mb_per_cpu=4000
    envmodules:
        "GCC/13.3.0"
    shell:        
        """
        /user/brussel/109/vsc10945/home/scratch/Software/kmerGWAS/bin/build_kmers_table -l {input.kmers_list} -k 31 -a {input.combined_kmers} -o {params.outfile} &> {log}
        """

rule kmers_table_to_bed:
    input:
        kmers = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table.table"),
        kmers_list = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_to_combine")
    output:
        plink = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_to_plink.0.bed")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "kmers_tobed.log")
    params:
        phenotype = os.path.join(RESULTS_DIR, "kmer", "combined", "phenotype.pheno"),
        kmers_table = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table"),
        kmers_to_plink = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_to_plink")
    resources:
        cpus_per_task=4,
        partition="zen4,zen5_mpi",
        mem_mb_per_cpu=90000
    envmodules:
        "GCC/13.3.0"
    shell:
        """
        echo -e "accession_id\tphenotype_value" > {params.phenotype}
        awk '{{print $2 "\t" ($2~/M/?1:2)}}' {input.kmers_list} >> {params.phenotype}
        /user/brussel/109/vsc10945/home/scratch/Software/kmerGWAS/bin/kmers_table_to_bed -t {params.kmers_table} -k 31 -p {params.phenotype} --maf 0.05 --mac 3 -b 1000000000 -o {params.kmers_to_plink}  &> {log}
        """

rule plink_to_abyss:
    input:
        plink0 = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_to_plink.0.bed")
    output:
        assoc_file = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers.assoc"),
        sig_assoc_file = os.path.join(RESULTS_DIR, "kmer", "combined", "significant_kmers.assoc"),
        male_assoc_file = os.path.join(RESULTS_DIR, "kmer", "combined", "male_significant_kmers.assoc"),
        abyss_input = os.path.join(RESULTS_DIR, "kmer", "combined", "abyss.input")
    params:
        workdir = os.path.join(RESULTS_DIR, "kmer", "combined")
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=16000
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "plink.log") 
    envmodules:
        "PLINK/2.00a3.7-foss-2022a",
        "Python/3.10.4-GCCcore-11.3.0"
    shell:
        """
        cd {params.workdir}
        plink --noweb --bfile kmers_to_plink.0 --allow-no-sex --assoc --out kmers &> {log}
        awk '$9 < 0.005' kmers.assoc | sed -E 's/[[:space:]]+/\t/g' > {output.sig_assoc_file}
        awk '$10 == 0' kmers.assoc | sed -E 's/[[:space:]]+/\t/g' > {output.male_assoc_file}
        python3 /user/brussel/109/vsc10945/home/scratch/Snakemake/SexLinked_SNP/scripts/plink_to_abyss_kmers.py {output.male_assoc_file} {output.abyss_input}
        """

rule abyss:
    input: 
        abyss_input = os.path.join(RESULTS_DIR, "kmer", "combined", "abyss.input")
    output:  
        abyss_output = os.path.join(RESULTS_DIR, "kmer", "combined", "abyss.output")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "abyss.log")  
    envmodules:
        "ABySS/2.3.7-foss-2023a"
    shell:
        """
        ABYSS -k25 -c0 -e0 {input.abyss_input} -o {output.abyss_output}
        """

rule blast_db:
    input:
        genome = REFERENCE_GENOME
    output:
        os.path.join(REF_DIR, f"{os.path.basename(REFERENCE_GENOME)}.nin")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "blastdb.log")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        makeblastdb -in {input.genome} -dbtype nucl 2> {log}
        """

rule blast_kmers:
    input:
        abyss_output = os.path.join(RESULTS_DIR, "kmer", "combined", "abyss.output"),
        genome = REFERENCE_GENOME,
        blast_db = os.path.join(REF_DIR, f"{os.path.basename(REFERENCE_GENOME)}.nin")
    output:
        blasted = os.path.join(RESULTS_DIR, "kmer", "combined", "blast.out")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmers", "blastout.log") 
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    shell:
        """
        blastn -query {input.abyss_output} -db {input.genome} -outfmt 6 > {output.blasted}  2> {log}
        """