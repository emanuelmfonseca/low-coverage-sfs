#!/bin/bash

#SBATCH --job-name=dadi_inference_CEU_YRI_ANGSD
#SBATCH --output=outfiles/dadi_inference/%x-%A-%a.out
#SBATCH --error=outfiles/dadi_inference/%x-%A-%a.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --array=0-3

### no low cov correction for ANGSD
# micromamba activate low-coverage when submit job

depths=(3x 5x 10x 30x)

mkdir -p "inference_split_mig"
python 3_dadi_inference.py ${depths[${SLURM_ARRAY_TASK_ID}]}
