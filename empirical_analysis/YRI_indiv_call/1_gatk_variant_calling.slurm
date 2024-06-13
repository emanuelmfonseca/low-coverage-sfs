#!/bin/bash

#SBATCH --job-name=gatk_indiv_var_call_YRI_GRCh38_job_array
#SBATCH --output=outfiles/gatk_variant_call/%x-%A-%a.out
#SBATCH --error=outfiles/gatk_variant_call/%x-%A-%a.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --array=1-80
### 4 depths * 20 samples = 80 jobs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# list of sample ID for task array
filename=$(ls ../YRI/chr20_bam_files/*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)

output_prefix=$(basename "$filename") # remove the dir in path

# reference genome
ref_path="../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# make outdir name depending on input file
if [[ $output_prefix == *"10x"* ]]; then 
    gatk_outdir="gatk_called_vcf_10x"
elif [[ $output_prefix == *"5x"* ]]; then
    gatk_outdir="gatk_called_vcf_5x"
elif [[ $output_prefix == *"3x"* ]]; then
    gatk_outdir="gatk_called_vcf_3x"
else
    gatk_outdir="gatk_called_vcf_30x"
fi

# make dir for output files
mkdir -p $gatk_outdir

# gatk
module load gatk
gatk --java-options "-Xmx4G" HaplotypeCaller -R $ref_path -L chr20 --native-pair-hmm-threads 4 \
-I $filename -O "${gatk_outdir}/${output_prefix}_chr20.vcf.gz" -ERC GVCF
