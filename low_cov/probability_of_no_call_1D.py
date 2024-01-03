import numpy
from dadi.LowCoverage import partitions_and_probabilities

def probability_of_no_call_1D(coverage_distribution, n_sequenced, Fx):
    """
    Calculate the probability of no genotype call for all allele frequencies.
    
    Args:
        coverage_distribution (numpy.ndarray): Coverage distribution.
        n_sequenced (int): Number of sequenced haplotypes.
    
    Returns:
        numpy.ndarray: Array containing the probability of no genotype call for each allele frequency.
    """
    # Extract partitions and their probabilities for a given number of samples
    partitions, partitions_probabilities = partitions_and_probabilities(n_sequenced, 'genotype', Fx)
    
    # Array of depths corresponding to coverage_distribution
    depths = numpy.arange(len(coverage_distribution[0]))
    
    # Create an empty array to store the final probabilities of no genotype calling
    all_prob_nocall = numpy.empty(n_sequenced + 1)
        
    for allele_freq, (partitions_, part_probs) in enumerate(zip(partitions, partitions_probabilities)):
        prob_nocall = 0
        
        for part, part_prob in zip(partitions_, part_probs):
            # Number of homozygous ref, heterozygous, and homozygous alt
            num_heterozygous, num_hom_alt = part.count(1), part.count(2)
            
            # Probability of getting no reads for the homozygous alt
            P_case0 = (
                coverage_distribution[1][0]**num_hom_alt *
                numpy.sum(coverage_distribution[1] * 0.5**depths)**num_heterozygous
            )
            
            # Probability of getting one read for the homozygous alt
            # P_case1a: Probability of 1 read in homozygous alt and 0 in heterozygous
            P_case1a = (
                num_hom_alt * coverage_distribution[1][1] * coverage_distribution[1][0]**(num_hom_alt - 1) *
                numpy.sum(coverage_distribution[1] * 0.5**depths)**num_heterozygous
            )
            
            # P_case1b: Probability of 1 read in heterozygous and 0 in homozygous alt
            P_case1b = (
                coverage_distribution[1][0]**num_hom_alt *
                numpy.sum(coverage_distribution[1] * 0.5**depths)**(num_heterozygous - 1) *
                num_heterozygous * numpy.sum(depths * coverage_distribution[1] * 0.5**depths)
            )
            
            # Calculate the probability of no call
            prob_nocall += part_prob * (P_case0 + P_case1a + P_case1b)
        
        all_prob_nocall[allele_freq] = prob_nocall
    
    return all_prob_nocall