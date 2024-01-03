import numpy
from dadi.LowCoverage import partitions_and_probabilities, flatten_nested_list, simulate_reads, subsample_genotypes_1D

def simulate_genotype_calling(coverage_distribution, allele_frequency, n_sequenced, n_subsampling, number_simulations, Fx):
    """
    Simulate the calling algorithm for alleles of a given frequency.
    
    Args:
        coverage_distribution (list): Coverage distribution for each population.
        allele_frequency (list): True allele frequency in sequenced samples for each population.
        n_sequenced (list): Number of sequenced haplotypes for each population.
        n_subsampling (list): Number of haplotypes after subsampling to account for missing data for each population.
        number_simulations (int): Number of loci to simulate for.
    
    Returns:
        numpy.ndarray: Frequency spectrum for n_subsampling haplotypes with observed allele frequencies resulting from the calling process.
    """
    # Initialize an array to record the allele frequencies from each simulation
    output_freqs = numpy.zeros(([x + 1 for x in n_subsampling]))
    
    # Record the number of individuals in each population
    pop_n_sequenced = [x // 2 for x in n_sequenced]
    
    # Extract partitions and their probabilities for each allele frequency value and population
    population_partitions = [partitions_and_probabilities(n, 'allele_frequency', Fx_, af) for (af, n, Fx_) in zip(allele_frequency, n_sequenced, Fx)]
    
    # Flatten the partitions and their probabilities for all populations
    combined_partitions = flatten_nested_list([_[0] for _ in population_partitions], '+')
    combined_part_probabilities = flatten_nested_list([_[1] for _ in population_partitions], '*')
    
    # Iterate through partitions and their probabilities
    for partition, partition_probability in zip(combined_partitions, combined_part_probabilities):
        nsimulations = number_simulations * partition_probability
        
        # Generate reads for all loci for this aggregate partition
        n_ref, n_alt = simulate_reads(coverage_distribution, partition, pop_n_sequenced, int(nsimulations))
        
        # Keep only loci identified as polymorphic
        t_alt = numpy.sum(n_alt, axis=1)
        n_ref, n_alt = n_ref[t_alt >= 2], n_alt[t_alt >= 2]
        
        # Update the allele frequency spectrum
        output_freqs.flat[0] += numpy.sum(t_alt < 2)
        
        # Calculate genotype calls for remaining loci
        genotype_calls = numpy.empty(n_ref.shape, dtype=int)
        genotype_calls[(n_ref == 0) & (n_alt == 0)] = 99 
        genotype_calls[(n_ref > 0) & (n_alt == 0)] = 0
        genotype_calls[(n_ref > 0) & (n_alt > 0)] = 1
        genotype_calls[(n_ref == 0) & (n_alt > 0)] = 2
        
        # Split genotype calls by population
        splits = numpy.cumsum(pop_n_sequenced)[:-1]
        split_genotype_calls = numpy.split(genotype_calls, splits, axis=1)
        
        # Handle subsampling to account for missing data
        split_nind_called = [numpy.sum(genotype_calls != 99, axis=1) for genotype_calls in split_genotype_calls]
        split_enough_calls = [ind_called >= n_subsampling_ // 2 for (ind_called, n_subsampling_) in zip(split_nind_called, n_subsampling)]
        all_enough_calls = numpy.logical_and.reduce(split_enough_calls)
        split_genotype_calls = [genotype_calls[all_enough_calls] for genotype_calls in split_genotype_calls]
        
        # Record loci without enough calls
        output_freqs.flat[0] += numpy.sum(all_enough_calls == False)
        
        # Calculate called allele frequencies for each population
        called_freqs = numpy.empty((len(split_genotype_calls[0]), len(n_sequenced)), int)
        for pop_ii, (n_subsampling_ii, n_sequenced_ii, genotype_calls) in enumerate(zip(n_subsampling, n_sequenced, split_genotype_calls)):
            if n_subsampling_ii != n_sequenced_ii:
                genotype_calls = subsample_genotypes_1D(genotype_calls, n_subsampling_ii)
            called_freqs[:,pop_ii] = numpy.sum(genotype_calls, axis=1)
        
        # Use the histogramdd function to generate the frequency spectrum for these genotype calls
        binning = [numpy.arange(n_subsampling_ii + 2) - 0.5 for n_subsampling_ii in n_subsampling]
        called_fs, _ = numpy.histogramdd(called_freqs, bins=binning)
        
        # Update the output frequency spectrum
        output_freqs += called_fs
        
    return output_freqs / numpy.sum(output_freqs)
