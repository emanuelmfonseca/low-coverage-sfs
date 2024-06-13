dir="/xdisk/rgutenk/lnt/projects/low_coverage"
COUNTER=1
for pop in YRI CEU; do
    for depth in 3x 5x 10x 30x; do
        if [ ${pop} == "YRI" ]; then
            echo "$(ls ../YRI/chr20_bam_files/*${depth}*.bam | head -n 10)" > pop_${COUNTER}_${depth}_bam.filelist
        else
            echo "$(ls ../CEU_YRI/chr20_bam_files/*${depth}*.bam | head -n 10)" > pop_${COUNTER}_${depth}_bam.filelist
        fi
    done
    COUNTER=$[$COUNTER +1]
done
