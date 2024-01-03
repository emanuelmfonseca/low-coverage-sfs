import numpy
import dadi
import warnings
from dadi.LowCoverage import part_inbreeding_probability

def partitions_and_probabilities(n_sequenced, partition_type, Fx=0, allele_frequency=None):
    """
    Generate allele count partitions and calculate their probabilities for a given number of alleles.
    
    Args:
        n_sequenced (int): The number of sequenced haplotypes.
        partition_type (str): Type of partition to generate ("allele_frequency" or "genotype").
        allele_frequency (int): The total number of alleles (required only for "allele_frequency" partition_type).
    
    Returns:
        tuple: A tuple containing two elements:
            - list: A list of all possible allele count partitions.
            - numpy.ndarray: An array containing the probability of each partition.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
         
        if partition_type == 'allele_frequency':
            # Check if n_sequenced is even, as partitions are calculated for even numbers
            if n_sequenced % 2 != 0:
                raise ValueError("Genotype partitions can only be calculated for an even number of haplotypes")
            
            # Generate partitions
            partitions = dadi.Numerics.cached_part(allele_frequency, n_sequenced/2)
            
            if Fx == 0:
                # Calculate partition probabilities
                partition_ways = numpy.array([numpy.exp(dadi.Numerics.multinomln([part.count(0), part.count(1), part.count(2)])) * 2 ** part.count(1) for part in partitions])
                partition_probabilities = partition_ways / numpy.sum(partition_ways)

            else:
                partition_probabilities = part_inbreeding_probability(partitions, Fx)
        
        elif partition_type == 'genotype':
            # Generate an array of allele counts from 0 to n_sequenced
            allele_counts = numpy.arange(n_sequenced + 1)
            
            # Generate partitions
            partitions = [dadi.Numerics.cached_part(allele_count, n_sequenced / 2) for allele_count in allele_counts]
            
            if Fx == 0:
                # Calculate partition probabilities using multinomial likelihood
                partition_ways = [
                    [
                        numpy.exp(dadi.Numerics.multinomln([part.count(0), part.count(1), part.count(2)])) * 2 ** part.count(1)
                        for part in parts
                    ]
                    for parts in partitions
                ]
                
                # Normalize partition probabilities
                partition_ways_sum = [[numpy.sum(part)] if len(part) > 1 else part for part in partition_ways]
                partition_probabilities = [numpy.array(pw) / numpy.array(pwb) for pw, pwb in zip(partition_ways, partition_ways_sum)]

            else:
                    partition_probabilities = [part_inbreeding_probability(part, Fx) for part in partitions]
         
        else:
            raise ValueError("Invalid partition_type. Use 'allele_frequency' or 'genotype'.")
        
        return partitions, partition_probabilities