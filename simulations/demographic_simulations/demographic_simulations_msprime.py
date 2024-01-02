# Import necessary libraries
import msprime
import os
import math
from pathlib import Path
import random

# Function for msprime simulation
def msprime_simulation(working_directory, model, params, values, nref, ns, ploidy, mu, seq_length, n_snps, recomb, rep, coverage, coverage_distribution, seed):
    
    # Extract model parameters
    model_params = {param: float(values[index]) for index, param in enumerate(params)}
    
    # Define demography based on the specified model
    if model == 'exp_growth_two_epochs':
        time_pop_change = 2 * nref * model_params['T1']
        ncurr = nref * model_params['nu1']
        r_curr= math.log(ncurr/nref)/time_pop_change
        demography = msprime.Demography()
        demography.add_population(name="pop1", description="population 1", initial_size=ncurr, growth_rate=r_curr)
        demography.add_population_parameters_change(time=time_pop_change, initial_size=nref, growth_rate=0)
    elif model == 'split_mig_two_pops':
        ncurr_pop1 = nref * model_params['nu1']
        ncurr_pop2 = nref * model_params['nu2']
        time_split = 2 * nref * model_params['T1']
        demography = msprime.Demography()
        demography.add_population(name="pop1", description="population 1", initial_size=ncurr_pop1)
        demography.add_population(name="pop2", description="population 2", initial_size=ncurr_pop2)
        demography.add_population(name="anc", description="ancestral population", initial_size=nref)
        demography.add_population_split(time=time_split, derived=["pop1", "pop2"], ancestral="anc")
    
    # Create necessary directories
    if not os.path.isdir(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'Genomic_files', 'Replicate_' + str(rep), 'Reference_genome')):
        os.makedirs(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'Genomic_files', 'Replicate_' + str(rep), 'Reference_genome'))
    
    if not os.path.isdir(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime')):
        os.makedirs(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime'))
    
    for replicate in range(1, rep + 1):
        if not os.path.isdir(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'Genomic_files', 'Replicate_' + str(replicate))):
            os.makedirs(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'Genomic_files', 'Replicate_' + str(replicate)))
    
    # Create an empty SNP file
    with open(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'Genomic_files', 'Replicate_' + str(rep), 'Reference_genome', \
    'SNPs_reference_genome_replicate' + str(rep) + '.txt'), 'w') as output:
        pass
        

    # Generate msprime tree sequences
    i = 0
    s = 0 
    while i < n_snps:
        # Choose the appropriate model and simulate ancestry
        if model == 'exp_growth_two_epochs' or model == 'inst_bottleneck_two_epochs':
            ts = msprime.sim_ancestry(samples={"pop1": ns[0]}, ploidy=ploidy, demography=demography, 
                                    sequence_length=1000, recombination_rate=recomb, random_seed=seed+i)
        elif model == 'split_mig_two_pops':
            ts = msprime.sim_ancestry(samples={"pop1": ns[0], "pop2": ns[1]}, ploidy=ploidy, demography=demography, 
                                    sequence_length=1000, recombination_rate=recomb, random_seed=seed+i)
    
        # Add mutations
        mts = msprime.sim_mutations(ts, rate=mu, discrete_genome=False) 
        
        # Save one SNP per simulation
        with open(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime', 'simulation_msprime_' + model.lower() + '_replicate' + str(rep) + '_prov.vcf'), "w") as vcf_file:
            mts.write_vcf(vcf_file, contig_id="0")
        
        # Handle the initial iteration
        if i == 0:
            fid = Path(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime', 'simulation_msprime_' + model.lower() + '_replicate' + str(rep) + '.vcf'))
            fid = open(fid, 'w')
            fid.close()
            
            with open(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime', 'simulation_msprime_' + model.lower() + '_replicate' + str(rep) + '_prov.vcf'), "r") as vcf_file:
                lines = vcf_file.readlines()
                for line in lines:
                    if line[0] == '#':
                        fid = open(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime', 'simulation_msprime_' + model.lower() + '_replicate' + str(rep) + '.vcf'), 'a')
                        silent = fid.write(line)
                        fid.close()
            
        save_snp = True
        with open(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime', 'simulation_msprime_' + model.lower() + '_replicate' + str(rep) + '_prov.vcf'), "r") as vcf_file:
            lines = vcf_file.readlines()
            random.seed(i)
            random.shuffle(lines)
            for line in lines:
                if line[0] != '#' and save_snp:
                    s += 500
                    l = line.split('\t')
                    l[1] = str(s)
                    fid = open(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime', 'simulation_msprime_' + model.lower() + '_replicate' + str(rep) + '.vcf'), 'a')
                    silent = fid.write('\t'.join(l))
                    fid.close()
                    i += 1
                    save_snp = False
            
        ii = 0
        with open(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'Genomic_files', 'Replicate_' + str(rep), 'Reference_genome', 'SNPs_reference_genome_replicate' + str(rep) + '.txt'), 'a') as output:
            for var in mts.variants():
                if ii == 0:
                    alleles = var.genotypes
                    alleles = [str(i) for i in alleles]
                    v = [str(s), '\t'.join(alleles)]
                    v = [str(i) for i in v]
                    output.write('\t'.join(v) + '\n')
                    ii += 1

    # Remove temporary file
    os.remove(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime', 'simulation_msprime_' + model.lower() + '_replicate' + str(rep) + '_prov.vcf'))
    
    nind = 0
    with open(os.path.join(working_directory, coverage_distribution, 'coverage_' + str(coverage), 'VCF_msprime', 'popfile.txt'), "w") as popfile:
        for pop, i in enumerate(ns):
            for ii in range(i):
                popfile.write('tsk_' + str(nind) + ' ' + 'pop' + str(pop + 1) + '\n')
                nind += 1
