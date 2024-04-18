#!/bin/bash

#SBATCH --job-name=angsd_YRI
#SBATCH --output=outfiles/angsd/%x-%A-%a.out
#SBATCH --error=outfiles/angsd/%x-%A-%a.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=05:00:00
#SBATCH --array=0-3

# program executable path
angsd='/xdisk/rgutenk/lnt/projects/low_coverage/software/angsd/angsd'
realSFS='/xdisk/rgutenk/lnt/projects/low_coverage/software/angsd/misc/realSFS'

# ancestral ref path
ancestral='../reference/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_20.fa'

# make result path
mkdir -p angsd_results

depths=(3x 5x 10x 30x)
depth=${depths[${SLURM_ARRAY_TASK_ID}]}

# run ANGSD
$angsd -bam ${depth}_bam.filelist -doSaf 1 -anc $ancestral -r chr20: -GL 2 -P 4 -minMapQ 1 -minQ 30 -out "angsd_results/${depth}"
$realSFS "angsd_results/${depth}.saf.idx" -maxIter 100 -P 4 > "angsd_results/${depth}.sfs.mle.txt"
