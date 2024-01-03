import numpy
from dadi.LowCoverage import partitions_and_probabilities

def calling_error_matrix(coverage_distribution, n_subsampling, Fx=0):
    """
        Calculate the calling error matrix based on the coverage distribution and subsampling.
        
        Args:
            coverage_distribution (numpy.ndarray): Coverage distribution.
            n_subsampling (int): Number of haplotypes to subsample.
        
        Returns:
            numpy.ndarray: Calling error matrix.
        """
    # Extract partitions and their probabilities for a given number of samples
    partitions, partitions_probabilities = partitions_and_probabilities(n_subsampling, 'genotype', Fx)
    
    # Array of depths corresopnding to coverage distribution
    depths = coverage_distribution[0][1:]
    
    # Probability that a heterozygote is called incorrectly
    # Note that zero reads would be no call, so it isn't included here
    coverage_distribution_ = [x/coverage_distribution[1][1:].sum() for x in coverage_distribution[1][1:]]
    prob_het_err = numpy.sum(coverage_distribution_ * 0.5**depths)
    
    # Transformation matrix
    trans_matrix = numpy.zeros((n_subsampling+1,n_subsampling+1))
    for allele_freq, (partitions_, part_probs) in enumerate(zip(partitions, partitions_probabilities)):
        for part, part_prob in zip(partitions_, part_probs):
            # Extract the number of heterozygous
            n_heterozygous = part.count(1)
            
            # For each possible number of heterozygous errors
            for n_error in range(n_heterozygous+1):
                # Probability that n_error heterozygous are called incorrectly out of heterozygous
                p_nerr = ssd.binom.pmf(n_error, n_heterozygous, prob_het_err)
                
                # Potential numbers of references erros, alternative erros, and the net change in allele frequency
                n_ref = numpy.arange(n_error+1)
                n_alt = n_error - n_ref
                net_change = n_alt - n_ref
                
                # Probabilities of each possible number of reference erros
                p_nref = ssd.binom.pmf(n_ref, n_error, 0.5)
                
                # Allele frequencies after errors
                afs_after_error = allele_freq + net_change
                
                # Record where error alleles end up
                trans_matrix[allele_freq, afs_after_error] += part_prob * p_nerr * p_nref
    
    return trans_matrix
