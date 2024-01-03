import numpy
import math
from scipy.special import comb

def probability_enough_individuals_covered(coverage_distribution, n_sequenced, n_subsampling):
    """
    Calculate the probability of having enough individuals covered to obtain n_subsampling successful genotypes.
    
    Args:
        coverage_distribution (numpy.ndarray): Coverage distribution.
        n_sequenced (int): Number of sequenced haplotypes.
        n_subsampling (int): Number of haplotypes to subsample.
    
    Returns:
        float: Probability of having enough individuals covered.
      """
    # Initialize the probability of having enough individuals covered
    prob_enough_individuals_covered = 0
        
    # Use math.ceil to round up, as you need a minimum of n_subsampling//2 covered individuals
    for covered in range(int(math.ceil(n_subsampling/2)), n_sequenced//2+1):
        # Calculate the probability of having enough individuals covered
        prob_enough_individuals_covered += (
            coverage_distribution[1][0]**(n_sequenced//2 - covered) *
            numpy.sum(coverage_distribution[1][1:]) ** covered *
            comb(n_sequenced//2, covered)
        )
    
    return prob_enough_individuals_covered