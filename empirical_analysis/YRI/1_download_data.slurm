#!/bin/bash

#SBATCH --job-name=download_YRI_30x_GRCh38
#SBATCH --output=outfiles/%x-%j.out
#SBATCH --error=outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lnt@arizona.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00

mkdir cram_files
cd cram_files

# get just the NA... sample name from popfile.txt
# cat popfile.txt | cut -f1 > file.txt
# grep -f file.txt all_YRI_aligned_wgs_GRCh38.tsv | grep cram | cut -f1 | sort > download_paths.txt

# download all files listed in the text file
while read line; do wget $line -nv; done < ../download_paths.txt