# command to make popfile.txt

# get dir where the BAM file for each pop is from
YRI_dir="../YRI/chr20_bam_files"
CEU_dir=chr20_bam_files"

# grab NA????? sample name
ls ${YRI_dir}/NA*3x*.bam | head -n 10 | grep -o 'NA[0-9]\+' | uniq > popfile.txt
# add pop name after a tab character
sed -i "s/$/\tYRI/" popfile.txt 

# repeat for CEU and append to popfile.txt
ls ${CEU_dir}/NA*3x*.bam | head -n 10 | grep -o 'NA[0-9]\+' | uniq >> popfile.txt
sed -i "s/$/\tCEU/" popfile.txt

# remove the extra CEU column in the YRI lines
# either cut for awk works
# has to save to a temp file first to not overwriting popfile.txt 
# (which is in the buffer anyway and won't work)
cut -f1,2 popfile.txt > temp.txt && mv temp.txt popfile.txt
# awk '{print $1,$2}' popfile.txt > temp.txt && mv temp.txt popfile.txt