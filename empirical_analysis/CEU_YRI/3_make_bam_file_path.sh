YRI_dir="../YRI/chr20_bam_files"
CEU_dir="chr20_bam_files"

for depth in 3x 5x 10x 30x; do
    echo "$(ls ${YRI_dir}/NA*${depth}*.bam | head -n 10)" > bam_file_path_${depth}.txt
    echo "$(ls ${CEU_dir}/NA*${depth}*.bam)" >> bam_file_path_${depth}.txt
done
