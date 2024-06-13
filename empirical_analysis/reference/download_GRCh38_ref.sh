# download source for reference genome Human GRCh38
# ref: https://github.com/dellytools/delly/issues/192
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.dict
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai

# download source for the corresponding ancestor sequence
# https://useast.ensembl.org/info/genome/compara/ancestral_sequences.html
wget https://ftp.ensembl.org/pub/release-110/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz