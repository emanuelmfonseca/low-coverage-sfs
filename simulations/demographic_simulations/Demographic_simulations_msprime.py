# Import necessary libraries
import os
import sys
import shutil

from Demographic_simulations_msprime import msprime_simulation

# Set the working directory for demographic models
working_directory = '/Users/emanuelmfonseca/projects/LowCoverageGenomes/Analyses/Demographic_models'

# Set parameters for simulations
replicates = 25
params = ['nu1', 'nu2', 'T1']
nu1 = [10, 1]
nu2 = ['NA', 1]
T1 = [0.1, 0.1]
SFS_dim = ['1D', '2D']

# Define coverage-related parameters
coverage_distributions = ['heterogeneous_coverage']
coverages = [3, 5, 10, 30]

# Iterate through specified parameters
for index, pop_size in enumerate(eval(params[0])):
    
    # Determine the demographic model and sample size based on SFS dimension
    if pop_size > 1 and SFS_dim[index] == '1D':
        model = 'exp_growth_two_epochs'
        ns = [20]
    elif SFS_dim[index] == '2D':
        model = 'split_mig_two_pops'
        ns = [10, 10]
    
    # Prepare parameters for simulation
    params_ = []
    values_ = []
    for param in params:
        s = str(eval(param)[index])
        if not s == 'NA':
            params_.append(param)
            values_.append(s)
    
    # Create a suffix for the model directory
    suffix = [params_[index2] + '_' + values_[index2] for index2, k in enumerate(params_)]
    suffix = '_'.join(suffix)
    
    # Define the model directory path
    model_directory = os.path.join(working_directory, SFS_dim[index] + '_' + model + '_' + suffix )
        
    # Iterate through coverage distributions
    for coverage_distribution in coverage_distributions:
        coverage = coverages[0]
        
        # Iterate through replicates
        for replicate in range(1, replicates + 1):
            
            # Create directories if they do not exist
            if not os.path.isdir(model_directory):
                os.makedirs(model_directory)
            
            if not os.path.isdir(os.path.join(model_directory, coverage_distribution, 'coverage_' + str(coverage), 'Genomic_files', 'Replicate_' + str(replicate))):
                os.makedirs(os.path.join(model_directory, coverage_distribution, 'coverage_' + str(coverage), 'Genomic_files', 'Replicate_' + str(replicate)))
        
            # Perform demographic simulation
            msprime_simulation(working_directory=model_directory,
                                model=model,
                                params=params_,
                                values=values_,
                                nref=1e4,
                                ns=ns,
                                ploidy=2,
                                mu=1e-8,
                                n_snps=20000,
                                recomb=1e-8,
                                rep=replicate,
                                coverage=coverage,
                                coverage_distribution=coverage_distribution,
                                seed=replicate)
    
    # Copy genomic files for different coverages
    coverages_ = coverages[1:]
    src_dir = os.path.join(model_directory, coverage_distribution, 'coverage_' + str(coverage))
    for coverage_ in coverages_:
        dst_dir = os.path.join(model_directory, coverage_distribution, 'coverage_' + str(coverage_))
        shutil.copytree(src_dir, dst_dir)
