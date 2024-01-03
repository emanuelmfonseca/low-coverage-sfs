# Import necessary libraries
from dadi.LowCoverage.LowCoverage import *
import dadi
import nlopt
import glob
import os
import sys
from pathlib import Path
import timeout_decorator

# Define a demographic model for bottleneck growth with inbreeding
def bottlegrowth_inbreeding(params, ns, pts):
    """
    Instantaneous size change followed by exponential growth.

    params = (nuB, nuF, T, F)
    ns = (n1,)

    nu1: Ratio of population size after instantaneous change to ancient
         population size
    nu2: Ratio of contemporary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations)
    n1: Number of samples in the resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1, nu2, T, F = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    nu_func = lambda t: nu1 * numpy.exp(numpy.log(nu2/nu1) * t/T)
    phi = dadi.Integration.one_pop(phi, xx, T, nu_func)
    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))

    return fs

# Timeout-decorated optimization function
@timeout_decorator.timeout(300)  # Set the timeout value in seconds
def optimize_with_timeout(p0, data_fs, demo_model_ex, grids, lower_bounds, upper_bounds):
    popt, ll_model = dadi.Inference.opt(p0, data_fs, demo_model_ex, grids,
                                        lower_bound=lower_bounds,
                                        upper_bound=upper_bounds,
                                        algorithm=nlopt.LN_COBYLA,
                                        maxeval=1000, verbose=0)
    return popt, ll_model

# Command line arguments
models = [sys.argv[1]]
coverage = sys.argv[2]
mode = sys.argv[3]
dataset = sys.argv[4]
nseq = [int(x) for x in sys.argv[5].split(',')]
nsub = [int(x) for x in sys.argv[6].split(',')]
coverage_distribution = sys.argv[7]

# Loop through models
for model in models:

    # Set working directory
    working_directory = '/xdisk/rgutenk/emanuelmfonseca/projects/LowCoverageGenomes/Analyses'

    # Determine simulation directory based on dataset
    if dataset == 'angsd':
        simulation_directory = os.path.join(working_directory, 'Demographic_models', model, coverage_distribution,
                                            'coverage_' + str(coverage), f'SFS_files_{dataset}')
    elif dataset == 'msprime' or dataset == 'SLiM':
        simulation_directory = os.path.join(working_directory, 'Demographic_models', model, coverage_distribution,
                                            'coverage_' + str(coverage), f'VCF_{dataset}')
    else:
        simulation_directory = os.path.join(working_directory, 'Demographic_models', model, coverage_distribution,
                                            'coverage_' + str(coverage), f'VCF_files_{dataset}')

    # Change to simulation directory
    os.chdir(simulation_directory)

    # Handle different file extensions based on dataset
    if dataset == 'angsd':
        files = glob.glob('*.fs')
    else:
        files = glob.glob('*.vcf')

    # Create lists for population IDs
    pop_ids = []
    all_pop_ids = []

    # Process population file for non-'angsd' datasets
    if dataset != 'angsd':
        popfile = os.path.join(simulation_directory, 'popfile.txt')
        with open(popfile) as popfile_infile:
            for line in popfile_infile:
                pop = line.split(' ')[1].rstrip()
                all_pop_ids.append(pop)

        # Deduplicate population IDs
        for pop_id in all_pop_ids:
            if pop_id not in pop_ids:
                pop_ids.append(pop_id)

    # Define optimization and individual directories
    optimization_directory = os.path.join(working_directory, 'Demographic_inference', model, coverage_distribution,
                                          'coverage_' + str(coverage), 'Optimized_parameters')
    individual_directory = os.path.join(working_directory, 'Demographic_inference', model, coverage_distribution,
                                        'coverage_' + str(coverage), 'Individual_parameters')

    # Create directories if they do not exist
    if not os.path.isdir(optimization_directory):
        os.makedirs(optimization_directory, exist_ok=True)

    if not os.path.isdir(individual_directory):
        os.makedirs(individual_directory, exist_ok=True)

    # Define headers based on the demographic model
    if 'two_epochs' in model:
        header1 = 'Model\tCoverage\tReplicate\tll_model\tnu1\tT\ttheta\n'
    elif 'split' in model:
        header1 = 'Model\tCoverage\tReplicate\tll_model\tnu1\tnu2\tT\tm\ttheta\n'
    elif 'bottlegrowth' in model:
        header1 = 'Model\tCoverage\tReplicate\tll_model\tnu1\tnu2\tT\tF\ttheta\n'

    # Define the best_params_file
    best_params_file = f'Demographic_params_{model}_coverage_{coverage}_{coverage_distribution}_{mode}_{dataset}_snps.txt'

    # Write header to the file
    fid = Path(os.path.join(optimization_directory, best_params_file))
    fid = open(fid, 'w')
    silent = fid.write(header1)
    fid.close()

    # Loop through files in the simulation directory
    index_rep = 0
    for file in files:
        # Split information from the file name
        info = file.split('_')

        index_rep += 1

        # Determine the replicate based on dataset
        if dataset == 'msprime' or dataset == 'SLiM' or dataset == 'angsd':
            rep = info[-1].split('.')[0]
            replicate = int(''.join(filter(str.isdigit, rep)))
        else:
            replicate = info[2]

        # Load data spectrum based on dataset
        if dataset == 'angsd':
            data_fs = dadi.Spectrum.from_file(os.path.join(simulation_directory, file))
        else:
            ss = {'pop1': int(nsub[0] // 2)}
            if 'split' in model:
                ss['pop2'] = int(nsub[1] // 2)

            # Create data dictionary from VCF file
            data_dict = dadi.Misc.make_data_dict_vcf(os.path.join(simulation_directory, file), popfile, subsample=ss)

            # Add additional information to the data dictionary
            for chrom_pos in data_dict:
                data_dict[chrom_pos]['outgroup_allele'] = data_dict[chrom_pos]['segregating'][0]
                data_dict[chrom_pos]['outgroup_context'] = data_dict[chrom_pos]['segregating'][0]

            # Create data spectrum from the data dictionary
            data_fs = dadi.Spectrum.from_data_dict(data_dict, pop_ids, nsub)

        # Define grids for optimization
        grids = [max(nsub) + 20, max(nsub) + 30, max(nsub) + 40]

        # Define parameters, bounds, and initial values based on the demographic model
        if 'two_epochs' in model:
            demo_model_ex = dadi.Numerics.make_extrap_func(dadi.Demographics1D.growth)
            param_names = ['nu1', 'T']
            lower_bounds = [1e-05, 1e-05]
            upper_bounds = [50, 1]
            params = [1, 0.05]
            header2 = '# Log(likelihood)\tnu1\tT\ttheta\n'
        elif 'split' in model:
            demo_model_ex = dadi.Numerics.make_extrap_func(dadi.Demographics2D.split_mig)
            param_names = ['nu1', 'nu2', 'T', 'm']
            lower_bounds = [1e-05, 1e-05, 1e-05, 0]
            upper_bounds = [50, 50, 1, 0]
            params = [1, 1, 0.05, 0]
            header2 = '# Log(likelihood)\tnu1\tnu2\tT\tm\ttheta\n'
        elif 'bottlegrowth' in model:
            demo_model_ex = dadi.Numerics.make_extrap_func(bottlegrowth_inbreeding)
            param_names = ['nu1', 'nu2', 'T', 'F']
            lower_bounds = [1e-05, 1e-05, 1e-05, 1e-5]
            upper_bounds = [50, 50, 1, 0.99999]
            params = [1, 1, 0.05, 0.05]
            header2 = '# Log(likelihood)\tnu1\tnu2\tT\tF\ttheta\n'

        # Define optimization parameters file
        opt_params_file = f'Demographic_params_{model}_coverage_{coverage}_{coverage_distribution}_{mode}_{dataset}_Replicate_{replicate}_Demographic_params.InferDM.txt'

        # Write header to the file
        fid = Path(os.path.join(individual_directory, opt_params_file))
        fid = open(fid, 'w')
        silent = fid.write(header2)
        fid.close()

        # Initialize variables for optimization
        res = []
        rep = 0
        nrep = 100

        # Loop through replicates for optimization
        while rep < nrep:
            print(f'Model: {model}; Coverage: {coverage}; Iteration: {rep}; Replicate {index_rep}')

            # Perturb parameters and perform optimization (with timeout)
            p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                                          lower_bound=lower_bounds)
            try:
                popt, ll_model = optimize_with_timeout(p0, data_fs, demo_model_ex, grids, lower_bounds, upper_bounds)
            except timeout_decorator.timeout_decorator.TimeoutError:
                print("Optimization timed out.")
                continue

            # Generate model spectrum and calculate optimal sfs scaling
            model_fs = demo_model_ex(popt, nsub, grids)
            theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, data_fs)

            # Write results to individual parameters file
            fid = open(os.path.join(individual_directory, opt_params_file), 'a')
            res = [ll_model] + list(popt) + [theta0]
            res = [str(i) for i in res]
            silent = fid.write('\t'.join(res) + '\n')
            fid.close()

            rep += 1

        # Retrieve optimization results
        llik = []
        opts = []
        with open(os.path.join(individual_directory, opt_params_file)) as input_file:
            infile = input_file.readlines()
            for line in infile:
                if line[0] == '-':
                    parameters = line.rstrip().split('\t')
                    llik.append(float(parameters[0]))
                    opts.append(parameters)

        # Identify the best parameters
        best_param_index = llik.index(max(llik))
        bestfit_params = opts[best_param_index]

        # Prepare results for writing to the main parameters file
        best_fit_res = [str(model), str(coverage), str(replicate)]
        silent = [best_fit_res.append(str(param)) for param in bestfit_params]

        # Write results to the main parameters file
        fid = open(os.path.join(optimization_directory, best_params_file), 'a')
        silent = fid.write('\t'.join(best_fit_res) + '\n')
        fid.close()
