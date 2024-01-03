# Function to generate simulated genomes
generate_genomes <- function(working_directory, model, nsnps, coverage, rep, coverage_distribution, seed) {
    
    bp <- c("C", "G", "A", "T")  # Nucleotide base pairs
    
    # Read simulations
    simulation <- read.table(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "Reference_genome", paste0("SNPs_reference_genome_replicate", rep, ".txt")))
    simulation <- simulation[1:nsnps, ]  # Extract N SNPs
    position <- round(simulation[, 1])  # Extract position information
    
    # Create diploid individuals
    alleles <- data.frame(t(simulation[, c(2:ncol(simulation))]))
    allele1_i <- seq(1, nrow(alleles), 2)
    allele2_i <- seq(2, nrow(alleles), 2)
    
    # Number of individuals
    num_inds <- nrow(alleles)/2
    
    # Simulate a random genome
    set.seed(seed * 64)
    reference_genome <- sample(bp, (max(position) + 500), replace = TRUE)
    genome_size <- length(reference_genome)  # Genome size
    reference_genome[position] <- rep("A", nsnps)
    reference_genome <- paste0(reference_genome, collapse = "")
    reference_genome <- substring(reference_genome, seq(1, genome_size, 70), seq(70, genome_size, 70))
    
    # Save the genome as a fasta file
    write(">Reference", file = file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "Reference_genome", paste0("reference_genome_replicate", rep, ".fasta")), append = FALSE)
    
    # Save a non-sequential fasta file
    for (line in reference_genome) {
        write(line, file = file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "Reference_genome", paste0("reference_genome_replicate", rep, ".fasta")), append = TRUE)
    }
    
    # Generate genomes from simulations
    for (ind in 1:num_inds) {
        
        allele1 <- alleles[allele1_i[ind], ]
        allele2 <- alleles[allele2_i[ind], ]
        genotype <- allele1 + allele2
        set.seed(seed * 64)
        genome1 <- sample(bp, (max(position) + 500), replace = TRUE, prob = c(0.2, 0.2, 0.3, 0.3))
        genome2 <- genome1
        
        snps_seq1 <- character()
        snps_seq2 <- character()
        
        # Generate sequences based on genotypes
        for (snp in genotype) {
            if (snp == 0) {
                snps_seq1 <- c(snps_seq1, "A")
                snps_seq2 <- c(snps_seq2, "A")
            } else if (snp == 1) {
                bases <- sample(c("A", "C"))
                snps_seq1 <- c(snps_seq1, bases[1])
                snps_seq2 <- c(snps_seq2, bases[2])
            } else {
                snps_seq1 <- c(snps_seq1, "C")
                snps_seq2 <- c(snps_seq2, "C")
            }
        }
        
        # Apply generated sequences to genomes
        genome1[position] <- snps_seq1
        genome1 <- paste0(genome1, collapse = "")
        genome1 <- substring(genome1, seq(1, genome_size, 70), seq(70, genome_size, 70))
        
        genome2[position] <- snps_seq2
        genome2 <- paste0(genome2, collapse = "")
        genome2 <- substring(genome2, seq(1, genome_size, 70), seq(70, genome_size, 70))
        
        folder_to_save <- file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), paste0("Genome_", ind))
        
        # Create folder to save the simulated genome
        if (!dir.exists(folder_to_save)) {
            dir.create(folder_to_save)
        }
        
        # Save genomes as fasta files
        write(paste0(">Genome1_", ind), file = file.path(folder_to_save, paste0("genome_", ind, ".fasta")), append = FALSE)
        for (line in genome1) {
            write(line, file = file.path(folder_to_save, paste0("genome_", ind, ".fasta")), append = TRUE)
        }
        
        write(paste0(">Genome2_", ind), file = file.path(folder_to_save, paste0("genome_", ind, ".fasta")), append = TRUE)
        for (line in genome2) {
            write(line, file = file.path(folder_to_save, paste0("genome_", ind, ".fasta")), append = TRUE)
        }
    }
}
