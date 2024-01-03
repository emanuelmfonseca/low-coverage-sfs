# Set the working directory and path to save Slurm scripts
working_directory <- "/Users/emanuelmfonseca/projects/LowCoverageGenomes/Analyses/Demographic_models"
path_to_save_scripts <- "/Users/emanuelmfonseca/projects/LowCoverageGenomes/Scripts/DM_Slurm_scripts"

# List of models, coverages, modes, datasets, and coverage distributions
models <- list.files(working_directory, recursive = FALSE)
coverages <- c(3,5,10,30)
modes = c('distortion', 'no_distortion')
datasets = c('gatk', 'angsd')
coverage_distributions <- 'heterogeneous_coverage'

# Function to generate Slurm script
Slurm_script <- function(model, coverage, mode, dataset, nseq, nsub, coverage_distribution){
    
    # Construct the Slurm script content
    Demographic_script <- paste0("#!/bin/bash
#SBATCH --job-name=Job_DM_simulation_", paste(c(model, mode, dataset), collapse="_"), "
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

module load anaconda/2020
conda init bash
source ~/.bashrc

python /Users/emanuelmfonseca/projects/LowCoverageGenomes/Scripts/Analyses_scripts/demographic_inference.py ", paste(model, coverage, mode, dataset, nseq, nsub, coverage_distribution))
    
    return(Demographic_script)
}

# Loop through different combinations of parameters to generate Slurm scripts
for (model in models){
    for (coverage in coverages){
        for (mode in modes){
            for (dataset in datasets){
                for (coverage_distribution in coverage_distributions){
                    if (mode == 'distortion' && dataset == 'angsd'){
                        # Skip this combination
                        next
                    } else {
                        # Set nseq and nsub based on conditions
                        if (dataset == 'angsd' && unlist(strsplit(model, "_"))[1] == '1D'){
                            nseq <- 40
                            nsub <- 40
                        } else if (dataset == 'angsd' && unlist(strsplit(model, "_"))[1] == '2D'){
                            nseq <- c(20, 20)
                            nsub <- c(20, 20)
                        } else if (dataset != 'angsd' && unlist(strsplit(model, "_"))[1] == '1D') {
                            nseq <- 40
                            nsub <- 32
                        } else if (dataset != 'angsd' && unlist(strsplit(model, "_"))[1] == '2D'){
                            nseq <- c(20,20)
                            nsub <- c(16,16)
                        }
                        
                        # Generate Slurm script content and save to file
                        cat("", file=file.path(path_to_save_scripts, paste0(paste(c(model, 'coverage', coverage, mode, dataset, 'nseq', paste(nseq, collapse="_"), 'nsub', paste(nsub, collapse="_"), coverage_distribution), collapse="_"),'.slurm')), append=FALSE)
                        cat(Slurm_script(model, coverage, mode, dataset, paste(nseq, collapse=","), paste(nsub, collapse=","), coverage_distribution), file=file.path(path_to_save_scripts, paste0(paste(c(model, 'coverage', coverage, mode, dataset, 'nseq', paste(nseq, collapse="_"), 'nsub', paste(nsub, collapse="_"), coverage_distribution), collapse="_"),'.slurm')), append=TRUE)
                    }
                }
            }
        }
    }
}
