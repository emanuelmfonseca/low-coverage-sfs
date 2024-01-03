import numpy
import scipy.stats as ss
from dadi.LowCoverage import split_list_by_lengths

def simulate_reads(coverage_distribution, flattened_partition, pop_n_sequenced, number_simulations):
    """
    Simulate reads across a given genotype partition.
    
    Args:
        coverage_distribution (list): Coverage distribution for each population.
        flattened_partition (list): Flattened genotype partition.
        pop_n_sequenced (list): Number of sequenced haplotypes for each population.
        number_simulations (int): Number of simulations to perform.
    
    Returns:
        tuple: A tuple containing two arrays:
            - numpy.ndarray: Arrays of reference allele counts for each simulated individual.
            - numpy.ndarray: Arrays of alternative allele counts for each simulated individual.
    """
    flattened_partition = numpy.array(flattened_partition)
    
    # Split partition for different populations
    partition = split_list_by_lengths(flattened_partition, pop_n_sequenced)
    
    # Empty array for storing the simulated data, initialized with zeros
    simulated_coverage = numpy.zeros((number_simulations, len(flattened_partition)), dtype=int)
    
    # Create population breaks
    splits = numpy.concatenate((numpy.array([0]), numpy.cumsum(pop_n_sequenced)), axis=0)
    
    # Simulate reads for each population
    pops = list(coverage_distribution.keys())
    for i, _ in enumerate(partition):
        cov_distribution = coverage_distribution[pops[i]][1]
        cov_sampling = ss.rv_discrete(values=[numpy.arange(len(cov_distribution)), cov_distribution])
        coverages = cov_sampling.rvs(size=(number_simulations, len(partition[i])))
        simulated_coverage[:, splits[i]:splits[i+1]] = coverages
    
    # Initialize arrays for reference and alternative allele counts
    n_ref, n_alt = numpy.zeros((number_simulations, len(flattened_partition)), dtype=int), numpy.zeros((number_simulations, len(flattened_partition)), dtype=int)
    
    # Calculate reference and alternative allele counts based on genotype partition
    n_ref[:, flattened_partition == 0] = simulated_coverage[:, flattened_partition == 0]
    n_alt[:, flattened_partition == 2] = simulated_coverage[:, flattened_partition == 2]
    n_alt[:, flattened_partition == 1] = ss.binom.rvs(simulated_coverage[:, flattened_partition == 1], 0.5)
    n_ref[:, flattened_partition == 1] = simulated_coverage[:, flattened_partition == 1] - n_alt[:, flattened_partition == 1]
    
    return n_ref, n_alt