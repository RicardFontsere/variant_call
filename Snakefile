import os
import glob

# Load configuration
configfile: "master/config/config.yaml"

# Include rule modules
include: "rules/common.smk"
include: "rules/preprocessing.smk"
include: "rules/mapping.smk"
include: "rules/coverage.smk"
include: "rules/variants.smk"
include: "rules/gwas.smk"
include: "rules/kmers.smk"
include: "rules/kmer_gwas.smk"

# Main rule
rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "trimmed", "{sample}_{read}_trimmed.fq.gz"), sample=SAMPLES, read=["1", "2"]),
        expand(os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "mapped", "{sample}.sorted.bam.csi"), sample=SAMPLES),
        REF_INDEX_FILES,
        GENOME_LENGTH,
        WINDOWS_FILE,
        expand(os.path.join(RESULTS_DIR, "coverage", "{sample}_multicov.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "coverage", "malenormalized", "{sample}_multicov_normalized.txt"), sample=ML_SAMPLES),
        expand(os.path.join(RESULTS_DIR, "coverage", "femalenormalized", "{sample}_multicov_normalized.txt"), sample=FL_SAMPLES),
        os.path.join(RESULTS_DIR, "coverage", "coverage_males.txt"),
        os.path.join(RESULTS_DIR, "coverage", "coverage_females.txt"),
        os.path.join(RESULTS_DIR, "coverage", "male2femalecov.txt"),
        expand(os.path.join(RESULTS_DIR, "variants", "variants_{species}.vcf.gz"), species=SPECIES),
        expand(os.path.join(RESULTS_DIR, "variants", "variants_{species}.vcf.gz.tbi"), species=SPECIES),
        expand(os.path.join(RESULTS_DIR, "variants", "variants_filtered_{species}.vcf.gz.tbi"), species=SPECIES),
        expand(os.path.join(RESULTS_DIR, "variants", "samples", "{sample}.vcf.gz"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "variants", "samples", "{sample}_50kb.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "variants", "samples", "malenormalized", "{sample}_50kb_normalized.txt"), sample=ML_SAMPLES),
        expand(os.path.join(RESULTS_DIR, "variants", "samples", "femalenormalized", "{sample}_50kb_normalized.txt"), sample=FL_SAMPLES),
        os.path.join(RESULTS_DIR, "variants", "snpdensity_males.txt"),
        os.path.join(RESULTS_DIR, "variants", "snpdensity_females.txt"),
        os.path.join(RESULTS_DIR, "variants", "snpdensity_fc.txt"),
        # GWAS  
        expand(os.path.join(RESULTS_DIR, "GWAS", "GWAS_{species}.ped"), species=SPECIES),
        expand(os.path.join(RESULTS_DIR, "GWAS", "GWAS_{species}.map"), species=SPECIES),
        expand(os.path.join(RESULTS_DIR, "GWAS", "GWAS_plink_{species}.bed"), species=SPECIES),
        expand(os.path.join(RESULTS_DIR, "GWAS", "chr_map_{species}.txt"), species=SPECIES),
        expand(os.path.join(RESULTS_DIR, "GWAS", "GWAS_gemma_{species}.assoc.txt"), species=SPECIES),
        # K-mers
        expand(os.path.join(RESULTS_DIR, "kmer", "{sample}", "output_kmc_canon.kmc_pre"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "kmer", "{sample}", "output_kmc_all.kmc_pre"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "kmer", "{sample}", "kmers_with_strand"), sample=SAMPLES),
        os.path.join(RESULTS_DIR, "kmer", "combined", "combined_kmers"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_to_combine"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table.table"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table.names"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_to_plink.0.bed"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "kmers.assoc"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "abyss.input"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "abyss.output"),
        os.path.join(REF_DIR, f"{os.path.basename(REFERENCE_GENOME)}.nin"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "blast.out"),
        os.path.join(RESULTS_DIR, "kmer", "combined", "significant_kmers.assoc")
        ## K-mer GWAS
        #os.path.join(RESULTS_DIR, "kmerGWAS", "sign_kmers.assoc"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "most_assoc_kmers.list"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "most_assoc_kmers.table"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "kmers.fasta"),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "mdiscards", "{sample}_discard_1.fq.gz"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "mdiscards", "{sample}_discard_2.fq.gz"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_matched_1.fq"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_matched_2.fq"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_stats"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "fdiscards", "{sample}_discard_1.fq.gz"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "fdiscards", "{sample}_discard_2.fq.gz"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_matched_1.fq"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_matched_2.fq"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_stats"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "mmegahit", "{sample}", "{sample}.contigs.fa"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "fmegahit", "{sample}", "{sample}.contigs.fa"), sample=FL_SAMPLES),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "female_list.txt"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "male_list.txt"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "FilteredFemaleKmers.txt"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "FilteredMaleKmers.txt"),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mdiscards", "{sample}_discard_1.fq.gz"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mdiscards", "{sample}_discard_2.fq.gz"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mbbduk", "{sample}_matched_1.fq"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mbbduk", "{sample}_matched_2.fq"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mbbduk", "{sample}_stats"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mmegahit", "{sample}", "{sample}.contigs.fa"), sample=ML_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fdiscards", "{sample}_discard_1.fq.gz"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fdiscards", "{sample}_discard_2.fq.gz"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fbbduk", "{sample}_matched_1.fq"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fbbduk", "{sample}_matched_2.fq"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fbbduk", "{sample}_stats"), sample=FL_SAMPLES),
        #expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fmegahit", "{sample}", "{sample}.contigs.fa"), sample=FL_SAMPLES),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "all_female_contigs.fasta"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "reduced_female_contigs.fasta"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "all_male_contigs.fasta"),
        #os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "reduced_male_contigs.fasta")