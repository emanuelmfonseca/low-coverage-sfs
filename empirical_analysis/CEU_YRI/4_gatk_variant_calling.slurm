#!/bin/bash

#SBATCH --job-name=gatk_var_calling_CEU_YRI_GRCh38_job_array
#SBATCH --output=outfiles/gatk_variant_call/%x-%A-%a.out
#SBATCH --error=outfiles/gatk_variant_call/%x-%A-%a.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --array=0-3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load gatk
module load samtools
module load openjdk/19.0.1

date
# picard path
picard_path="/xdisk/rgutenk/lnt/projects/low_coverage/software/picard.jar"

i=${SLURM_ARRAY_TASK_ID}
depths=("3x" "5x" "10x" "30x")
depth=${depths[${i}]}

merge_dir='merged_bam_files'
mkdir -p $merge_dir

dedup_dir='dedup_bam_files'
mkdir -p $dedup_dir

gatk_outdir="gatk_called_vcf"
mkdir -p $gatk_outdir

# reference genome
ref_path="../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"

output_prefix=CEU_YRI_${depth}_GRCh38_chr20

# merge bam files
samtools merge -b bam_file_path_${depth}.txt ${merge_dir}/${output_prefix}_merged.bam

# mark duplicates
java -jar $picard_path MarkDuplicates \
            -I ${merge_dir}/${output_prefix}_merged.bam \
            -O ${dedup_dir}/${output_prefix}_merged_dedup.bam \
            -M ${dedup_dir}/${output_prefix}_merged_dedup_metrics.txt
                
# index
samtools index ${dedup_dir}/${output_prefix}_merged_dedup.bam


# gatk
# variant calling
gatk --java-options "-Xmx4G" HaplotypeCaller \
        -R $ref_path -L chr20 \
        --native-pair-hmm-threads $SLURM_CPUS_PER_TASK \
        -I ${dedup_dir}/${output_prefix}_merged_dedup.bam \
        -O ${gatk_outdir}/${output_prefix}_gatk_raw_variants.vcf
        
# select only SNP
gatk --java-options "-Xmx4G" SelectVariants \
        -R $ref_path \
        -V ${gatk_outdir}/${output_prefix}_gatk_raw_variants.vcf \
        --select-type-to-include SNP \
        -O ${gatk_outdir}/${output_prefix}_gatk_only_SNPs.vcf
        
# filter snps
gatk --java-options "-Xmx4G" VariantFiltration \
        -R $ref_path \
        -V "${gatk_outdir}/${output_prefix}_gatk_only_SNPs.vcf" \
        --filter-name 'QD_filter' --filter-expression 'QD<2.0' \
        --filter-name 'FS_filter' --filter-expression 'FS>60.0' \
        --filter-name 'MQ_filter' --filter-expression 'MQ<40.0' \
        --filter-name 'SOR_filter' --filter-expression 'SOR>10.0' \
        -O "${gatk_outdir}/${output_prefix}_gatk_filtered.vcf"

date
