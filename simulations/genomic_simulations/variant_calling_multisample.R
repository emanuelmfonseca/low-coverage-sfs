variant_calling_multisample <- function(working_directory, model, ns, coverage, rep, coverage_distribution){
    
    # Set paths to necessary tools
    angsd_path <- '[path_to_angsd]'
    gatk_path <- '[path_to_gatk]'
    bwa_path <- '[path_to_bwa]'
    samtools_path <- '[path_to_samtools]'
    picard_path <- '[path_to_picard]'
    
    # Calculate the number of individuals based on files in the working directory
    ninds <- length(list.files(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep)), pattern="Genome", full.names = TRUE))
    dimensions <- unlist(strsplit(model, '_'))[1]
    
    # Define key paths for reference and merged files
    ind_ref <- file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "Reference_genome")
    merged_files_path <- file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "Genomic_files", paste0("Replicate_", rep), "merged_files")
    
    # Delete .dict file if it exists
    if (file.exists(file.path(ind_ref, paste0("reference_genome_replicate", rep,".dict")))){
        file.remove(file.path(ind_ref, paste0("reference_genome_replicate", rep,".dict")))
    }
    
    # Prepare reference dictionary, fasta index, and bam index
    system(paste0("java -jar $PICARDJARPATH/", picard_path, " CreateSequenceDictionary R=", file.path(ind_ref, paste0("reference_genome_replicate", rep,".fasta")), " O=", file.path(ind_ref, paste0("reference_genome_replicate", rep,".dict"))))
    
    system(paste(samtools_path, "faidx", file.path(ind_ref, paste0("reference_genome_replicate", rep,".fasta"))))
    
    # Mark duplicates
    system(paste0("java -jar $PICARDJARPATH/", picard_path ," MarkDuplicates INPUT=", file.path(merged_files_path, 'merged.bam'), " OUTPUT=", file.path(merged_files_path, "dedup_aligned_reads.bam"), " METRICS_FILE=", file.path(merged_files_path, "metrics.txt")))
    
    system(paste(samtools_path, " index", file.path(merged_files_path, "dedup_aligned_reads.bam")))
    
    # Variant calling - GATK
    system(paste0(gatk_path, " --java-options '-Xmx4g' HaplotypeCaller -R ", file.path(ind_ref, paste0("reference_genome_replicate", rep,".fasta")), " -I ", file.path(merged_files_path, "dedup_aligned_reads.bam"), " -O ", file.path(merged_files_path, "gatk_raw_variants.vcf")))
    
    # Extract only SNPs
    system(paste0(gatk_path, " --java-options '-Xmx4g' SelectVariants -R ", file.path(ind_ref, paste0("reference_genome_replicate", rep,".fasta")), " -V ", file.path(merged_files_path, "gatk_raw_variants.vcf")," --select-type-to-include SNP -O ", file.path(merged_files_path, "gatk_only_SNPs_variants.vcf")))
    
    # Create folder to store the final vcf
    if (!dir.exists(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_gatk"))){
        dir.create(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_gatk"))
    }
    
    # Filter SNPs
    system(paste0(gatk_path, " --java-options '-Xmx4g' VariantFiltration -R ", file.path(ind_ref, paste0("reference_genome_replicate", rep,".fasta")), " -V ", file.path(merged_files_path, "gatk_only_SNPs_variants.vcf"), " --filter-name 'QD_filter' --filter-expression 'QD<2.0' --filter-name 'FS_filter' --filter-expression 'FS>60.0' --filter-name 'MQ_filter' --filter-expression 'MQ<40.0' --filter-name 'SOR_filter' --filter-expression 'SOR>10.0' -O ", file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_gatk", paste0("gatk_Replicate_", rep ,"_filtered.vcf"))))
    
    # Remove index files if they exist
    if (length(list.files(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_gatk"), pattern="*.idx")) > 0){
        file.remove(list.files(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_gatk"), pattern="*.idx", full.names = TRUE))
    }
    
    # Create a file to store population information
    d <- as.numeric(strsplit(dimensions, '')[[1]][1])
    cat("", file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_gatk", 'popfile.txt'), append=FALSE)
    cat("", file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_bcftools", 'popfile.txt'), append=FALSE)
    cat("", file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_freebayes", 'popfile.txt'), append=FALSE)
    ninds_ <- 1
    for (x in 1:d){
        for (i in 1:(ninds/d)){
            info <- paste(c('tsk_', ninds_,  " pop", x, "\n"), collapse="")
            cat(info, file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_gatk", 'popfile.txt'), append=TRUE)
            cat(info, file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_bcftools", 'popfile.txt'), append=TRUE)
            cat(info, file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage),  "VCF_files_freebayes", 'popfile.txt'), append=TRUE)
            ninds_ <- ninds_ + 1
        }
    }
    
    
    ## Estimating a SFS using ANGSD ##
    
    # Create folder to store final SFS
    if (!dir.exists(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd"))){
        dir.create(file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd"))
    }
    
    # Extract number of populations
    SFS_dim <- unlist(strsplit(model, "_"))[1]
    
    # Index the BAM files
    system(paste("for i in", file.path(merged_files_path, '*_sorted_aligned_reads.bam; do samtools index $i; done')))
    
    # Generate file lists
    if (SFS_dim == '1D'){
        system(paste("ls", file.path(merged_files_path, '*_sorted_aligned_reads.bam'), '>', file.path(merged_files_path, 'pop1.filelist')))
        
        # Filter out low-quality bases and reads
        system(paste(file.path(angsd_path, 'angsd'), '-bam', file.path(merged_files_path, 'pop1.filelist'), '-doSaf 1 -out', file.path(merged_files_path, paste0('angsd_SFS_replicate_', rep)), '-anc', file.path(ind_ref, paste0('reference_genome_replicate', rep,'.fasta')), '-GL 2 -P 4 -minMapQ 1 -minQ 30'))
        
        # Obtain maximum likelihood estimate of SFS using EM algorithm
        system(paste(file.path(angsd_path, 'misc', 'realSFS'), file.path(merged_files_path,paste0('angsd_SFS_replicate_',rep, '.saf.idx')), '-maxIter 100 -P 4 >', file.path(merged_files_path, paste0('angsd_SFS_final_replicate_',rep, '.txt'))))
        
        sfs <- read.table(file.path(merged_files_path, paste0("angsd_SFS_final_replicate_",rep, ".txt")), header=FALSE)
        
        cat("", file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd", paste0("angsd_SFS_replicate_", rep, ".fs")), append=FALSE)
        
        nind <- length(list.files(merged_files_path, pattern='_sorted_aligned_reads.bam$'))
        cat(paste((nind * 2 + 1), 'unfolded "pop1"\n'), file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd", paste0("angsd_SFS_replicate_", rep, ".fs")), append=TRUE)
        cat(paste0(paste0(sfs, collapse=' '), "\n"), file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd", paste0("angsd_SFS_replicate_", rep, ".fs")), append=TRUE)
        
        mask <- paste0(c(1, rep(0, (length(sfs)-2)), 1), collapse=" ")
        cat(paste0(mask, "\n", collapse=" "), file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd", paste0("angsd_SFS_replicate_", rep, ".fs")), append=TRUE)
        
    } else if (SFS_dim == '2D') {
        pop_index <- split(1:ninds, rep(1:length(ns), c(ns)))
        for (i in 1:length(ns)){
            system(paste("ls", file.path(merged_files_path, paste0('*_{',min(unlist(pop_index[i])), '..', max(unlist(pop_index[i])), '}_sorted_aligned_reads.bam', '>', file.path(merged_files_path, paste0('pop', i, '.filelist'))))))
            
            # Filter out low-quality bases and reads
            system(paste(file.path(angsd_path, 'angsd'), '-bam', file.path(merged_files_path, paste0('pop', i, '.filelist')), '-doSaf 1 -out', file.path(merged_files_path, paste0('angsd_SFS_replicate_',rep, '_pop_', i)), '-anc', file.path(ind_ref, paste0('reference_genome_replicate',rep,'.fasta')), '-GL 2 -P 4 -minMapQ 1 -minQ 30'))
            
        }
        
        # Obtain maximum likelihood estimate of SFS using EM algorithm
        system(paste(file.path(angsd_path, 'misc', 'realSFS'), paste(file.path(merged_files_path,paste0('angsd_SFS_replicate_',rep, '_pop_', 1:length(ns), '.saf.idx')), collapse=' '), '-maxIter 100 -P 4 >', file.path(merged_files_path, paste0('angsd_SFS_final_replicate_',rep, '.txt'))))
        
        sfs <- read.table(file.path(merged_files_path, paste0("angsd_SFS_final_replicate_",rep, ".txt")), header=FALSE)
        
        cat("", file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd", paste0("angsd_SFS_replicate_", rep, ".fs")), append=FALSE)
        
        nind <- length(list.files(merged_files_path, pattern='_sorted_aligned_reads.bam$'))
        cat(paste((ns[1] * 2 + 1), (ns[2] * 2 + 1), 'unfolded "pop1" "pop2"\n'), file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd", paste0("angsd_SFS_replicate_", rep, ".fs")), append=TRUE)
        cat(paste0(paste0(sfs, collapse=' '), "\n"), file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd", paste0("angsd_SFS_replicate_", rep, ".fs")), append=TRUE)
        
        mask <- paste0(c(1, rep(0, (length(sfs)-2)), 1), collapse=" ")
        cat(paste0(mask, "\n", collapse=" "), file=file.path(working_directory, model, coverage_distribution, paste0("coverage_", coverage), "SFS_files_angsd", paste0("angsd_SFS_replicate_", rep, ".fs")), append=TRUE)
        
    }
    
}
