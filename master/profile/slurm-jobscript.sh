#!/bin/bash
#SBATCH --job-name=job_%x.%A
#SBATCH --output=logs/%x.%A.out
#SBATCH --error=logs/%x.%A.err
#SBATCH --ntasks={resources.tasks}
#SBATCH --cpus-per-task={resources.cpus_per_task}
#SBATCH --mem-per-cpu={resources.mem_mb_per_cpu}M
#SBATCH --time={resources.time}

# Load modules
module load snakemake/7.22.0-foss-2022a
module load Java/11.0.20
module load FastQC/0.11.9-Java-11
module load Trimmomatic/0.39-Java-11
module load Bowtie2/2.4.5-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

# Run the command
{exec_job}

