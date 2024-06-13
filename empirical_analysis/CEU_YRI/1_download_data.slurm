#!/bin/bash

#SBATCH --job-name=download_CEU_30x_GRCh38
#SBATCH --output=outfiles/%x-%j.out
#SBATCH --error=outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lnt@arizona.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=8:00:00

## runtime: 05:15:56
# ignore the first line with column label
# tail -n +2 CEU_30x_GRCh38.tsv 

# sort, take the first 20 and only get first columns for download link
# tail -n +2 CEU_30x_GRCh38.tsv | sort | head -20 | cut -f1
# can skip sort if desired

# generate download path from CEU_30x_GRCh38.tsv list downloaded from 1K Genome website
# filter CEU > 1000 Genomes 30x on GRCh38 > Alignment > PCR-free high coverage (179 samples)
# tail -n +2 CEU_30x_GRCh38.tsv | head -20 | cut -f1 | sort > download_paths.txt

# Create a directory named fastq_files and change to that directory
mkdir -p cram_files; cd cram_files

# download all files listed in the text file
while read line; do wget $line -nv; done < ../download_paths.txt