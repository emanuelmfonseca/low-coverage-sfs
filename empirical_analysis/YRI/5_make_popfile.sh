# make popfile.txt
cd chr20_bam_files/ # or any data dir
# grab NA????? sample name
ls | grep -o 'NA[0-9]\+' | uniq > ../popfile.txt
# add pop name after a tab character
sed -i "s/$/\tYRI/" ../popfile.txt 
