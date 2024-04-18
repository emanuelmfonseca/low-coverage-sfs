#!/bin/bash

#SBATCH --job-name=cram_to_bam_job_array
#SBATCH --output=outfiles/cram_to_bam_array/%x-%A_%a.out
#SBATCH --error=outfiles/cram_to_bam_array/%x-%A_%a.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=05:00:00
#SBATCH --array=1-10

module load samtools
module load openjdk/19.0.1

# make dir for storing the bam files
mkdir -p bam_files
mkdir -p chr20_bam_files

# picard path
picard_path="/xdisk/rgutenk/lnt/projects/low_coverage/software/picard.jar"

# list of sample ID for task array
filename=$(ls cram_files/*.cram | sed -n ${SLURM_ARRAY_TASK_ID}p)

# remove ".final.cram" and "cram_files" dir to get output prefix
output_prefix=$(basename "${filename%.final*}")

# convert cram to bam
# each convertion takes about 1hr (tested for 10 cpus)
samtools view -b -o "bam_files/${output_prefix}.bam" $filename

# produce the index file with picard for each bam file
java -jar $picard_path BuildBamIndex -I "bam_files/${output_prefix}.bam" -O "bam_files/${output_prefix}.bai"

# filter for chr20 bam only
samtools view -b "bam_files/${output_prefix}.bam" chr20 > chr20_bam_files/${output_prefix}.30x.chr20.bam

# produce the index file with picard for each chr20 bam file
java -jar $picard_path BuildBamIndex -I chr20_bam_files/${output_prefix}.30x.chr20.bam -O chr20_bam_files/${output_prefix}.30x.chr20.bai

# produce subsampled files
samtools view -s 0.333 -b chr20_bam_files/${output_prefix}.30x.chr20.bam > chr20_bam_files/${output_prefix}.10x.chr20.bam
samtools view -s 0.166 -b chr20_bam_files/${output_prefix}.30x.chr20.bam > chr20_bam_files/${output_prefix}.5x.chr20.bam
samtools view -s 0.1 -b chr20_bam_files/${output_prefix}.30x.chr20.bam > chr20_bam_files/${output_prefix}.3x.chr20.bam

# produce index file with picard
java -jar $picard_path BuildBamIndex -I chr20_bam_files/${output_prefix}.10x.chr20.bam -O chr20_bam_files/${output_prefix}.10x.chr20.bai
java -jar $picard_path BuildBamIndex -I chr20_bam_files/${output_prefix}.5x.chr20.bam -O chr20_bam_files/${output_prefix}.5x.chr20.bai
java -jar $picard_path BuildBamIndex -I chr20_bam_files/${output_prefix}.3x.chr20.bam -O  chr20_bam_files/${output_prefix}.3x.chr20.bai

date