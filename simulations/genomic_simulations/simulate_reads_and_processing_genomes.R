# Function to simulate processing of genomes using in silico sequencing tools
simulate_reads_and_processing_genomes <- function(working_directory, model, coverage, read_length, rep, coverage_distribution) {
    
    # List genome directories for each replicate
    dirs <- list.files(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep)), pattern="Genome", full.names = TRUE)
    
    for (dir in dirs){
        
        # Extract individual ID from the directory name
        chr <- unlist(strsplit(dir, '/'))
        chr <- tail(chr, n=1)
        ind <- as.numeric(gsub('.*_(\\d+)$','\\1',chr))
        
        # Create abundance file for in silico sequencing
        ref1 <- paste0(c(paste0("Genome1_", ind), 0.5), collapse="\t")
        ref2 <- paste0(c(paste0("Genome2_", ind), 0.5), collapse="\t")
        
        write(ref1, file = file.path(dir, paste0("abundance_genome_", ind, ".txt")), append = FALSE)
        write(ref2, file = file.path(dir, paste0("abundance_genome_", ind, ".txt")), append = TRUE)
        
        # Path to the reference genome
        reference_genome <- file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "Reference_genome", paste0("reference_genome_replicate", rep, ".fasta"))
        
        # Calculate genome size
        ref_gen <- readLines(reference_genome)
        genome_size <- 0
        for (x in 2:length(ref_gen)){
            genome_size <- genome_size + nchar(ref_gen[x])
        }
        
        # Total number of reads to be simulated, considering coverage and read length
        if (coverage_distribution == 'heterogeneous_coverage'){
            coverage_ <- 0
            seed <- rep
            while (coverage_ < 1){
                set.seed(seed)
                coverage_ <- rnorm(1,coverage,0.1*coverage)
                seed <- seed + 1
            }
        } else if (coverage_distribution == 'uniform_coverage') {
            coverage_ <- coverage
        }
        
        # Adjust coverage for specific model conditions
        if (grepl('2D_het', model)){
            if (ind > 10){
                coverage_ = coverage_ / 2
            }
        }
        
        n_reads <- round(coverage_ * genome_size / read_length)
        
        # Set paths for tools and files
        ind_genome_file <- file.path(dir, paste0("genome_", ind, ".fasta"))
        abundance_file <- file.path(dir, paste0("abundance_genome_", ind, ".txt"))
        output_file <- file.path(dir, paste0("genome_", ind, ".fasta"))
        insilicoseq = '/home/u14/emanuelmfonseca/miniconda3/envs/env/bin/iss'
        bwa_path <- 'bwa'
        samtools_path <- 'samtools'
        path_softwares <- paste0(unlist(strsplit(working_directory, "/"))[2:6], collapse="/")
        picard_path <- file.path(path_softwares, "Softwares", "picard.jar")
        
        # Generate reads using in silico sequencing
        system(paste0(insilicoseq, " generate --model hiseq --genomes ",ind_genome_file, " --abundance_file ", abundance_file," --n_reads ", n_reads ," --cpus 4 --output ", output_file))
        
        # Generate index for the reference genome
        system(paste(bwa_path, "index", reference_genome))
        
        # Create folder to save bam files
        if (!dir.exists(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files"))){
            dir.create(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files"))
        }
        
        # Perform alignment
        system(paste0(bwa_path, " mem -M -R '@RG\\tID:tsk_",ind,"\\tLB:tsk_",ind,"\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:tsk_", ind, "' ", reference_genome, " ", file.path(dir, paste0("genome_", ind,".fasta_R1.fastq ")), file.path(dir, paste0("genome_", ind,".fasta_R2.fastq "))," > ", file.path(dir, paste0("genome_",ind,"_aligned_reads.sam"))))
        
        # Sort sam and convert to bam
        system(paste0("java -jar $PICARDJARPATH/", picard_path, " SortSam I=", file.path(dir, paste0("genome_",ind,"_aligned_reads.sam")), " O=", file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files", paste0("genome_",ind,"_sorted_aligned_reads.bam")), " SORT_ORDER=coordinate"))
        
    }
    
    # Collect information about the sorted alignment for each replicate
    sorted_alignment <- list.files(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files"), pattern="sorted_aligned")
    
    # Create read group file for merged bam file
    cat(NULL, file = file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files", 'rg.txt'), append = FALSE)
    
    # Add read group information to the read group file
    for (dir_ in dirs){
        i <- as.numeric(gsub('.*_(\\d+)$','\\1',dir_))
        rg <- paste0("@RG\tID:tsk_",i,"\tLB:tsk_",i,"\tPL:ILLUMINA\\tPM:HISEQ\tSM:tsk_", i,'\n')      
        cat(rg, file = file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files", 'rg.txt'), append = TRUE)
    }
    
    # Merge bam files for each replicate
    system(paste0(samtools_path, " merge -f ", file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files", 'merged.bam'), " ", paste0(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files", sorted_alignment), collapse=" ")))
    
}
