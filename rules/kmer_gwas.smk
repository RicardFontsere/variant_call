rule ref_kmer_list:
    input: 
        assoc_file = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers.assoc")
    output:
        sign_kmers = os.path.join(RESULTS_DIR, "kmerGWAS", "sign_kmers.assoc"),
        kmer_list = os.path.join(RESULTS_DIR, "kmerGWAS", "most_assoc_kmers.list")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "ref_kmer_list.log") 
    envmodules:
        "BBMap/39.01-GCC-11.3.0"
    shell:
        """
        mkdir -p $(dirname {output.sign_kmers})
        mkdir -p $(dirname {log})
        awk '$9 < 0.005' {input.assoc_file} | sed -E 's/[[:space:]]+/\t/g' > {output.sign_kmers}
        cat {output.sign_kmers} | cut -f 3 > {output.kmer_list}
        """

rule filter_kmers:
    input:
        kmer_list = os.path.join(RESULTS_DIR, "kmerGWAS", "most_assoc_kmers.list"),
        kmer_table_names = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table.names"),
        kmer_table_table = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table.table")
    output:
        assoc_table = os.path.join(RESULTS_DIR, "kmerGWAS", "most_assoc_kmers.table")
    params:
        kmer_table_base = os.path.join(RESULTS_DIR, "kmer", "combined", "kmers_table")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "filter_kmers.log") 
    shell:
        """
        /user/brussel/109/vsc10945/home/scratch/Software/kmerGWAS/bin/filter_kmers \
            -t {params.kmer_table_base} \
            -k {input.kmer_list} \
            -o {output.assoc_table}
        """

rule make_kmer_fasta:
    input:
        assoc_table = os.path.join(RESULTS_DIR, "kmerGWAS", "most_assoc_kmers.table")
    output:
        kmers_fasta = os.path.join(RESULTS_DIR, "kmerGWAS", "kmers.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "make_kmer_fasta.log")        
    shell:
        """
        tail -n +2 {input.assoc_table} \
        | cut -f 1 \
        | sed 's/^\\([^acgt]\\)/>kmer\\n\\1/' \
        | awk '/^>/{{sub(/>/,"&" ++i "_")}}1' > {output.kmers_fasta}
        """
    
rule male_bbduk:
    input:
        r1 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_1_trimmed.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_2_trimmed.fq.gz"),
        ref = os.path.join(RESULTS_DIR, "kmerGWAS", "kmers.fasta")
    output:
        mdiscard1 = os.path.join(RESULTS_DIR, "kmerGWAS", "mdiscards", "{sample}_discard_1.fq.gz"),
        mdiscard2 = os.path.join(RESULTS_DIR, "kmerGWAS", "mdiscards", "{sample}_discard_2.fq.gz"),
        mout1 = os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_matched_1.fq"),
        mout2 = os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_matched_2.fq"),
        mstats = os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_stats")
    params:
        k = 31
    envmodules:
        "BBMap/39.01-GCC-11.3.0"
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=40000
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "male_bbduk_{sample}.log") 
    shell:
        """
        mkdir -p {RESULTS_DIR}/kmerGWAS/mbbduk
        mkdir -p {RESULTS_DIR}/kmerGWAS/mdiscards
        mkdir -p $(dirname {log})
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.mdiscard1} \
            out2={output.mdiscard2} \
            outm={output.mout1} \
            outm2={output.mout2} \
            ref={input.ref} \
            k={params.k} \
            stats={output.mstats} \
            usejni=t
        """
rule female_bbduk:
    input:
        r1 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_1_trimmed.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_2_trimmed.fq.gz"),
        ref = os.path.join(RESULTS_DIR, "kmerGWAS", "kmers.fasta")
    output:
        fdiscard1 = os.path.join(RESULTS_DIR, "kmerGWAS", "fdiscards", "{sample}_discard_1.fq.gz"),
        fdiscard2 = os.path.join(RESULTS_DIR, "kmerGWAS", "fdiscards", "{sample}_discard_2.fq.gz"),
        fout1 = os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_matched_1.fq"),
        fout2 = os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_matched_2.fq"),
        fstats = os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_stats")
    params:
        k = 31
    envmodules:
        "BBMap/39.01-GCC-11.3.0"
    resources:
        cpus_per_task=4,
        mem_mb_per_cpu=40000
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "female_bbduk_{sample}.log") 
    shell:
        """
        mkdir -p {RESULTS_DIR}/kmerGWAS/fbbduk
        mkdir -p {RESULTS_DIR}/kmerGWAS/fdiscards
        mkdir -p $(dirname {log})
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.fdiscard1} \
            out2={output.fdiscard2} \
            outm={output.fout1} \
            outm2={output.fout2} \
            ref={input.ref} \
            k={params.k} \
            stats={output.fstats} \
            usejni=t    
        """

rule male_megahit:
    input:
        mout1 = os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_matched_1.fq"),
        mout2 = os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_matched_2.fq")
    output:
        mhitout = os.path.join(RESULTS_DIR, "kmerGWAS", "mmegahit", "{sample}", "{sample}.contigs.fa")
    params:
        outdir = os.path.join(RESULTS_DIR, "kmerGWAS", "mmegahit", "{sample}"),
        samples = "{sample}"
    resources:
        cpus_per_task=40,
        runtime=2040,
        mem_mb_per_cpu=4000 
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "{sample}_megahit.log") 
    shell:
        """
        /user/brussel/109/vsc10945/home/scratch/Software/MEGAHIT/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
            -1 {input.mout1} \
            -2 {input.mout2} \
            -f \
            -t 40 \
            -m 1 \
            --out-dir {params.outdir} \
            --out-prefix {params.samples} 2> {log}
        """

rule female_megahit:
    input:
        fout1 = os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_matched_1.fq"),
        fout2 = os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_matched_2.fq")
    output:
        fhitout = os.path.join(RESULTS_DIR, "kmerGWAS", "fmegahit", "{sample}", "{sample}.contigs.fa")
    params:
        outdir = os.path.join(RESULTS_DIR, "kmerGWAS", "fmegahit", "{sample}"),
        samples = "{sample}"
    resources:
        cpus_per_task=40,
        runtime=2040,
        mem_mb_per_cpu=4000 
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "{sample}_fmegahit.log") 
    shell:
        """
        /user/brussel/109/vsc10945/home/scratch/Software/MEGAHIT/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
            -1 {input.fout1} \
            -2 {input.fout2} \
            -f \
            -t 40 \
            -m 1 \
            --out-dir {params.outdir} \
            --out-prefix {params.samples} 2> {log}
        """

rule R_kmers:
    input:
        presab = os.path.join(RESULTS_DIR, "kmerGWAS", "most_assoc_kmers.table"),
        fstats = expand(os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk", "{sample}_stats"), sample=FL_SAMPLES),
        mstats = expand(os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk", "{sample}_stats"), sample=ML_SAMPLES),
    output:
        femaleIndstest = os.path.join(RESULTS_DIR, "kmerGWAS", "female_list.txt"),
        maleIndstest = os.path.join(RESULTS_DIR, "kmerGWAS", "male_list.txt"),
        filtered_female = os.path.join(RESULTS_DIR, "kmerGWAS", "FilteredFemaleKmers.txt"),
        filtered_male = os.path.join(RESULTS_DIR, "kmerGWAS", "FilteredMaleKmers.txt")
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "R_kmers.log")
    params:
        female_samples = " ".join(FL_SAMPLES),
        male_samples = " ".join(ML_SAMPLES),
        fbbduk = os.path.join(RESULTS_DIR, "kmerGWAS", "fbbduk"),
        mbbduk = os.path.join(RESULTS_DIR, "kmerGWAS", "mbbduk"),
        n_females = len(FL_SAMPLES),
        n_males = len(ML_SAMPLES)
    envmodules:
        "R/4.2.1-foss-2022a"
    resources:
        cpus_per_task=10,
        runtime=1200,
        mem_mb_per_cpu=4000
    shell:
        """
        Rscript --vanilla -e '
            # Read in data
            kmer_counts <- read.delim("{input.presab}", header=TRUE)
            # Create sample lists as data frames
            femaleInds <- data.frame(V1 = strsplit("{params.female_samples}", " ")[[1]])
            maleInds <- data.frame(V1 = strsplit("{params.male_samples}", " ")[[1]])
            # Write sample lists
            write.table(femaleInds, "{output.femaleIndstest}", row.names=F, col.names=F, quote=F)
            write.table(maleInds, "{output.maleIndstest}", row.names=F, col.names=F, quote=F)
            # Process female files
            female_files <- list.files(path = "{params.fbbduk}", pattern = "_stats")
            if(length(female_files) > 0) {{
                DF_Females <- read.delim(paste0("{params.fbbduk}/", female_files[1]), 
                                        header=FALSE, comment.char="#")
                colnames(DF_Females) <- c("KmerName", female_files[1], "perc")
                DF_Females$perc <- NULL
                if(length(female_files) > 1) {{
                    for (i in female_files[-1]) {{
                        df <- read.delim(paste0("{params.fbbduk}/", i), 
                                       header=FALSE, comment.char="#")
                        colnames(df) <- c("KmerName", i, "perc")
                        df$perc <- NULL
                        DF_Females <- merge(DF_Females, df, all=T, by="KmerName")
                    }}
                }}
                rownames(DF_Females) <- DF_Females$KmerName
                DF_Females$KmerName <- NULL
                DF_Females[is.na(DF_Females)] <- 0
                DF_Females$countszero <- rowSums(DF_Females == 0)
            }}
            # Process male files
            male_files <- list.files(path = "{params.mbbduk}", pattern = "_stats")
            if(length(male_files) > 0) {{
                DF_Males <- read.delim(paste0("{params.mbbduk}/", male_files[1]), 
                                      header=FALSE, comment.char="#")
                colnames(DF_Males) <- c("KmerName", male_files[1], "perc")
                DF_Males$perc <- NULL
                if(length(male_files) > 1) {{
                    for (i in male_files[-1]) {{
                        df <- read.delim(paste0("{params.mbbduk}/", i), 
                                       header=FALSE, comment.char="#")
                        colnames(df) <- c("KmerName", i, "perc")
                        df$perc <- NULL
                        DF_Males <- merge(DF_Males, df, all=T, by="KmerName")
                    }}
                }}
                rownames(DF_Males) <- DF_Males$KmerName
                DF_Males$KmerName <- NULL
                DF_Males[is.na(DF_Males)] <- 0
                DF_Males$countszero <- rowSums(DF_Males == 0)
            }}
            # Calculate total reads (adjust column indices based on number of samples)
            n_female_samples <- {params.n_females}
            n_male_samples <- {params.n_males}
            DF_Males$totalread <- rowSums(DF_Males[,1:n_male_samples])
            DF_Females$totalread <- rowSums(DF_Females[,1:n_female_samples])
            # Combine k-mer data
            DF_combKmer <- data.frame(DF_Females$totalread, DF_Females$countszero)
            DF_combKmer$kmer <- row.names(DF_Females)
            DF_temp <- data.frame(DF_Males$totalread, DF_Males$countszero)
            DF_temp$kmer <- row.names(DF_Males)
            DF_combKmer <- merge(DF_combKmer, DF_temp, all.x = T, all.y = T)
            # Handle NAs
            DF_combKmer$DF_Females.totalread[is.na(DF_combKmer$DF_Females.totalread)] <- 0
            DF_combKmer$DF_Males.totalread[is.na(DF_combKmer$DF_Males.totalread)] <- 0
            DF_combKmer$DF_Females.countszero[is.na(DF_combKmer$DF_Females.countszero)] <- n_female_samples
            DF_combKmer$DF_Males.countszero[is.na(DF_combKmer$DF_Males.countszero)] <- n_male_samples
            # Remove kmers present in all individuals
            DF_combKmer <- subset(DF_combKmer, 
                                 !(DF_combKmer$DF_Females.countszero == 0 & 
                                   DF_combKmer$DF_Males.countszero == 0))
            # Define sex-specific kmers
            # Female-specific: present in all females, absent in all males
            DF_fem_specific <- subset(DF_combKmer, 
                                     DF_combKmer$DF_Females.countszero == 0 & 
                                     DF_combKmer$DF_Males.countszero == n_male_samples)
            # Male-specific: present in all males, absent in all females
            DF_male_specific <- subset(DF_combKmer, 
                                      DF_combKmer$DF_Males.countszero == 0 & 
                                      DF_combKmer$DF_Females.countszero == n_female_samples)
            # Write output files
            write.table(DF_fem_specific, file = "{output.filtered_female}", 
                       sep = "\\t", quote = F, row.names = F)
            write.table(DF_male_specific, file = "{output.filtered_male}", 
                       sep = "\\t", quote = F, row.names = F)
        ' &> {log}
        """

rule male_bbduk_sexspecific:
    input:
        r1 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_1_trimmed.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_2_trimmed.fq.gz"),
        ref = os.path.join(RESULTS_DIR, "kmerGWAS", "kmers.fasta"),
        filtered_male = os.path.join(RESULTS_DIR, "kmerGWAS", "FilteredMaleKmers.txt")
    output: 
        mdiscard1 = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mdiscards", "{sample}_discard_1.fq.gz"),
        mdiscard2 = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mdiscards", "{sample}_discard_2.fq.gz"),
        mout1 = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mbbduk", "{sample}_matched_1.fq"),
        mout2 = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mbbduk", "{sample}_matched_2.fq"),
        mstats = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mbbduk", "{sample}_stats"),
        mhitout = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mmegahit", "{sample}", "{sample}.contigs.fa")
    params:
        filtered_male_kmers = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "{sample}_FilteredMaleKmers.fasta"),
        k = 31,
        outdir = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mmegahit", "{sample}"),
        samples = "{sample}"
    envmodules:
        "BBMap/39.01-GCC-11.3.0"
    resources:
        cpus_per_task=40,
        runtime=2040,
        mem_mb_per_cpu=4000 
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "male_bbdukss_{sample}.log")
    shell:
        """
        tail -n +2 {input.filtered_male} | cut -f1 | awk '{{print ">"$1}}' | grep -A 1 --no-group-separator -f - {input.ref} > {params.filtered_male_kmers} 2>> {log}
        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.mdiscard1} \
            out2={output.mdiscard2} \
            outm={output.mout1} \
            outm2={output.mout2} \
            ref={params.filtered_male_kmers} \
            k={params.k} \
            stats={output.mstats} \
            usejni=t >> {log} 2>&1
        /user/brussel/109/vsc10945/home/scratch/Software/MEGAHIT/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
            -1 {output.mout1} \
            -2 {output.mout2} \
            -f \
            -t 40 \
            -m 1 \
            --out-dir {params.outdir} \
            --out-prefix {params.samples} >> {log} 2>&1
        """

rule female_bbduk_sexspecific:
    input:
        r1 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_1_trimmed.fq.gz"),
        r2 = os.path.join(RESULTS_DIR, "trimmed", "{sample}_2_trimmed.fq.gz"),
        ref = os.path.join(RESULTS_DIR, "kmerGWAS", "kmers.fasta"),
        filtered_female = os.path.join(RESULTS_DIR, "kmerGWAS", "FilteredFemaleKmers.txt")
    output: 
        fdiscard1 = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fdiscards", "{sample}_discard_1.fq.gz"),
        fdiscard2 = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fdiscards", "{sample}_discard_2.fq.gz"),
        fout1 = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fbbduk", "{sample}_matched_1.fq"),
        fout2 = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fbbduk", "{sample}_matched_2.fq"),
        fstats = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fbbduk", "{sample}_stats"),
        fhitout = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fmegahit", "{sample}", "{sample}.contigs.fa")
    params:
        filtered_female_kmers = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "{sample}_FilteredFemaleKmers.fasta"),
        k = 31,
        outdir = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fmegahit", "{sample}"),
        samples = "{sample}"
    envmodules:
        "BBMap/39.01-GCC-11.3.0"
    resources:
        cpus_per_task=40,
        runtime=2040,
        mem_mb_per_cpu=4000 
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "female_bbdukss_{sample}.log") 
    shell:
        """
        tail -n +2 {input.filtered_female} | cut -f1 | awk '{{print ">"$1}}' | grep -A 1 --no-group-separator -f - {input.ref} > {params.filtered_female_kmers} 2>> {log}

        bbduk.sh \
            in={input.r1} \
            in2={input.r2} \
            out={output.fdiscard1} \
            out2={output.fdiscard2} \
            outm={output.fout1} \
            outm2={output.fout2} \
            ref={params.filtered_female_kmers} \
            k={params.k} \
            stats={output.fstats} \
            usejni=t  >> {log} 2>&1
        /user/brussel/109/vsc10945/home/scratch/Software/MEGAHIT/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
            -1 {output.fout1} \
            -2 {output.fout2} \
            -f \
            -t 40 \
            -m 1 \
            --out-dir {params.outdir} \
            --out-prefix {params.samples} >> {log} 2>&1
        """

rule reduce_male_contigs:
    input:
        mhitout = expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "mmegahit", "{sample}", "{sample}.contigs.fa"), sample=ML_SAMPLES)
    output:
        mcontigs = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "all_male_contigs.fasta"),
        mcontigsreduced = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "reduced_male_contigs.fasta")
    envmodules:
        "CD-HIT/4.8.1-GCC-13.3.0"
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=500,
        runtime=1200
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "reduced_male.log")
    shell:
        """
        cat {input.mhitout} > {output.mcontigs}
        cd-hit-est -i {output.mcontigs} -o {output.mcontigsreduced} -c 0.98 -n 10 -d 0 -M 16000 -T 40 >> {log} 2>&1
        """

rule reduce_female_contigs:
    input:
        fhitout = expand(os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "fmegahit", "{sample}", "{sample}.contigs.fa"), sample=FL_SAMPLES)
    output:
        fcontigs = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "all_female_contigs.fasta"),
        fcontigsreduced = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "reduced_female_contigs.fasta")
    envmodules:
        "CD-HIT/4.8.1-GCC-13.3.0"
    resources:
        cpus_per_task=20,
        mem_mb_per_cpu=500,
        runtime=1200
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "reduced_female.log")
    shell:
        """
        cat {input.fhitout} > {output.fcontigs}
        cd-hit-est -i {output.fcontigs} -o {output.fcontigsreduced} -c 0.98 -n 10 -d 0 -M 16000 -T 40 >> {log} 2>&1
        """ 
        
rule blast_male_contigs:
    input:
        mcontigsreduced = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "reduced_male_contigs.fasta"),
        genome = REFERENCE_GENOME,
        blast_db = os.path.join(REF_DIR, f"{os.path.basename(REFERENCE_GENOME)}.nin")
    output:
        blastout = os.path.join(RESULTS_DIR, "kmerGWAS", "blast", "males_kmers_blast.out")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    resources:
        cpus_per_task=40,
        mem_mb_per_cpu=16000,
        runtime=600
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "male_blast.log")
    shell:
        """
        blastn -query {input.mcontigsreduced} -db {input.genome} -outfmt 6 -out {output.blastout} -num_threads 39 >> {log} 2>&1
        """

rule blast_female_contigs:
    input:
        fcontigsreduced = os.path.join(RESULTS_DIR, "kmerGWAS", "sexspecific", "reduced_female_contigs.fasta"),
        genome = REFERENCE_GENOME,
        blast_db = os.path.join(REF_DIR, f"{os.path.basename(REFERENCE_GENOME)}.nin")
    output:
        blastout = os.path.join(RESULTS_DIR, "kmerGWAS", "blast", "females_kmers_blast.out")
    envmodules:
        "BLAST+/2.16.0-gompi-2024a"
    resources:
        cpus_per_task=40,
        mem_mb_per_cpu=16000,
        runtime=600
    log:
        os.path.join(RESULTS_DIR, "logs", "kmer_gwas", "female_blast.log")
    shell:
        """
        blastn -query {input.fcontigsreduced} -db {input.genome} -outfmt 6 -out {output.blastout} -num_threads 39 >> {log} 2>&1
        """ 