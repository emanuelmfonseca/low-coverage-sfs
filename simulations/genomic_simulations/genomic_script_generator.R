# Define working directory
working_directory <- "/Users/emanuelmfonseca/projects/LowCoverageGenomes/Analyses/Demographic_models"
path_to_save_scripts <- "/Users/emanuelmfonseca/projects/LowCoverageGenomes/Scripts/simulate_reads_and_processing_genomes"

# List of available models in the working directory
models_ <- list.files(working_directory, recursive = FALSE)
nsnps_ <- 20000
read_length_ <- 126
coverages_ <- c(3,5,10,30)
coverage_distributions_ <- 'heterogeneous_coverage'

# Function to generate an R script for simulation and variant calling
R_script <- function(model, nsnps, coverages, coverage_distributions, read_length){
  
  script <- paste0('library(parallel)

# Define working directory (duplicated - consider removing)
working_directory <- "/Users/emanuelmfonseca/projects/LowCoverageGenomes/Analyses/Demographic_models"

# Load necessary R scripts for simulation and variant calling
source("/Users/emanuelmfonseca/simulations/genomic_simulations/generate_genomes.R")
source("/Users/emanuelmfonseca/simulations/genomic_simulations/simulate_reads_and_processing_genomes.R")
source("/Users/emanuelmfonseca/simulations/genomic_simulations/variant_calling_multisample.R")

# Set parameters for simulation
model_ <- "', model, '"
nsnps_ <- ', nsnps, '
read_length_ <- ', read_length, '
coverage_ <- ', coverages, '
coverage_distribution_ <- "', coverage_distributions, '"

# Determine the number of dimensions for the site frequency spectrum (SFS)
SFS_dim <- unlist(strsplit(model_, "_"))[1]

# Set the number of samples based on SFS dimensions
if (SFS_dim == "1D"){
    ns_i = 20
} else if (SFS_dim == "2D"){
    ns_i = c(10, 10)
}

# Determine the number of replicates
nrep <- length(list.files(file.path(working_directory, model_, coverage_distribution_, paste0("coverage_", coverage_), "Genomic_files"), pattern="Replicate", recursive = FALSE))

# Function to generate whole genomes in parallel
generate_whole_genomes <- function(i, model_i, nsnps_i, coverage_i, coverage_distribution_i, read_length_i){
    
    # Generate individual genomes for each replicate
    generate_genomes(working_directory,
                     model = model_i,
                     nsnps = nsnps_i,
                     coverage = coverage_i,
                     rep = i,
                     coverage_distribution = coverage_distribution_i,
                     seed = 64)
    
    # Simulate reads from each individual and combine all of them in a bam file
    simulate_reads_and_processing_genomes(working_directory,
                                model = model_i,
                                coverage = coverage_i,
                                read_length = read_length_i,
                                rep = i,
                                coverage_distribution = coverage_distribution_i)
    
    # Variant calling using BWA and GATK / SFS using ANGSD
    variant_calling_multisample(working_directory,
                                model = model_i,
                                ns = ns_i,
                                coverage = coverage_i,
                                rep = i,
                                coverage_distribution = coverage_distribution_i)
}

# Run replicates in parallel using mclapply
mclapply(1:nrep, function(i, model_i, nsnps_i, coverage_i, coverage_distribution_i, read_length_i) generate_whole_genomes(i, model_,  nsnps_, coverage_, coverage_distribution_, read_length_),
         mc.cores=detectCores())

# Return the generated R script
')

  return(script)
}

# Function to generate a Slurm script for job submission
Slurm_script <- function(model, coverages, coverage_distributions){
  
  script <- paste0('#!/bin/bash
#SBATCH --job-name=Job_Sim_var_calling_', model,'
### SLURM reads %x as the job name and %j as the job ID
#SBATCH --output=%x-%j.out
#SBATCH --account=[user_defined_parameter]
#SBATCH --partition=[user_defined_parameter]
### REQUIRED. Set the number of cores that will be used for this job.
#SBATCH --ntasks=[user_defined_parameter]
### REQUIRED. Set the memory required for this job.
#SBATCH --mem=[user_defined_parameter]
### REQUIRED. Specify the time required for this job, hhh:mm:ss
#SBATCH --time=[user_defined_parameter]

# Load necessary modules
module load samtools/1.10
module load bwa/0.7.17
module load picard/2.23.4
module load gatk/4.2.5.0
module load bcftools/1.10.2
module load freebayes/1.3.2
module load R --no-save

# Run the generated R script for simulation and variant calling
R --no-save < /Users/emanuelmfonseca/projects/LowCoverageGenomes/Scripts/Simulate_genomes_and_var_calling/', paste(c(model, 'coverage', coverages, 'coverage_distribution', coverage_distributions),collapse='_'), '.R')

# Return the Slurm script
return(script)

}

# Loop through models, coverages, and coverage distributions to generate R and Slurm scripts
for (model in models_){
    for (coverage in coverages_){
        for (coverage_distribution in coverage_distributions_){
            # Generate and append R script to a file
            cat('', file=file.path(path_to_save_scripts, paste0(paste(c(model,'coverage', coverage, 'coverage_distribution', coverage_distribution), collapse='_'), '.R')), append=FALSE)
            cat(R_script(model, nsnps_, coverage, coverage_distribution, read_length_), file=file.path(path_to_save_scripts, paste0(paste(c(model,'coverage', coverage, 'coverage_distribution', coverage_distribution), collapse='_'), '.R')), append=TRUE)
            
            # Generate and append Slurm script to a file
            cat("", file=file.path(path_to_save_scripts, paste0(paste(c(model,'coverage', coverage, 'coverage_distribution', coverage_distribution), collapse='_'), '.slurm')), append=FALSE)
            cat(Slurm_script(model, coverage, coverage_distribution), file=file.path(path_to_save_scripts, paste0(paste(c(model,'coverage', coverage, 'coverage_distribution', coverage_distribution), collapse='_'), '.slurm')), append=TRUE)
        }
    }
}
