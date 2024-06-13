#!/bin/bash

#SBATCH --job-name=dadi_inference_CEU_YRI_gatk_chr20_job_array
#SBATCH --output=outfiles/dadi_inference_CEU_YRI/%x-%A-%a.out
#SBATCH --error=outfiles/dadi_inference_CEU_YRI/%x-%A-%a.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --array=0-3
#### for future run should run the low_cov case separately --> 4 instead of 2 jobs

# micromamba activate low-coverage when submit job

mkdir -p inference_split_mig
depth=("3x" "5x" "10x" "30x")
python dadi_inference.py ${depth[${SLURM_ARRAY_TASK_ID}]}
