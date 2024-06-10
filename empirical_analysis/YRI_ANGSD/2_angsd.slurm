#!/bin/bash

#SBATCH --job-name=angsd_YRI
#SBATCH --output=outfiles/angsd_GL_3/%x-%A-%a.out
#SBATCH --error=outfiles/angsd_GL_3/%x-%A-%a.err
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
ref_fa=/xdisk/rgutenk/lnt/projects/low-coverage-sfs/empirical_analysis/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

# make result path
mkdir -p angsd_results_SOAPsnp
mkdir -p angsd_results_SYK

depths=(3x 5x 10x 30x)
depth=${depths[${SLURM_ARRAY_TASK_ID}]}

# run ANGSD
# $angsd -bam ${depth}_bam.filelist -doSaf 1 -anc $ancestral -r chr20: -GL 2 -P 4 \
#        -minMapQ 1 -minQ 30 -out "angsd_results/${depth}"

# $angsd -bam ${depth}_bam.filelist -doSaf 1 -anc $ancestral -r chr20: -GL 1 -P 4 \
#         -minMapQ 30 -minQ 20 -minInd 16 -out "angsd_results/${depth}"

# $angsd -bam ${depth}_bam.filelist -doSaf 1 -anc $ancestral -r chr20: -GL 2 -P 10 \
#         -minMapQ 1 -minQ 30 -out "angsd_results/${depth}"


$angsd -bam ${depth}_bam.filelist -doSaf 1 -anc $ancestral -r chr20 -GL 3 -P 5 \
     -ref $ref_fa -out "angsd_results_SOAPsnp/${depth}"
     
$angsd -bam ${depth}_bam.filelist -GL 3 -out "angsd_results_SOAPsnp/${depth}" -doGlf 1 -r chr20
     
# $angsd -bam ${depth}_bam.filelist -doSaf 1 -anc $ancestral -r chr20 -GL 4 -P 5 \
#      -doCounts 1 -out "angsd_results_SYK/${depth}"


# $realSFS "angsd_results/${depth}.saf.idx" -P 10 > "angsd_results/${depth}.sfs.mle.txt"

$realSFS "angsd_results_SOAPsnp/${depth}.saf.idx" -P 5 -r chr20 > "angsd_results_SOAPsnp/${depth}.sfs.mle.txt"

# $realSFS "angsd_results_SYK/${depth}.saf.idx" -P 5 -r chr20 > "angsd_results_SYK/${depth}.sfs.mle.txt"
