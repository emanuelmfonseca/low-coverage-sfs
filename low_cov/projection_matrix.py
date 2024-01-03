import numpy
import dadi
from dadi.LowCoverage import partitions_and_probabilities, projection_inbreeding

def projection_matrix(n_sequenced, n_subsampling, F):
    """
    Create a projection matrix for down-sampling haplotypes.
    
    Args:
        n_sequenced (int): Number of haplotypes sequenced.
        n_subsampling (int): Number of haplotypes to project down to.
    
    Returns:
        numpy.ndarray: Projection matrix.
    
    """
    # Create an empty matrix to store the projection
    projection_matrix = numpy.empty((n_sequenced + 1, n_subsampling + 1))
    
    # Calculate the projection for each allele frequency
    for allele_freq in range(n_sequenced + 1):
        if F != 0:
            partitions, partition_probabilities = partitions_and_probabilities(n_sequenced, 'allele_frequency', F, allele_freq)
            
            proj = numpy.zeros_like(range(n_subsampling+1), dtype=float)
            for partition, part_prob in zip(partitions, partition_probabilities):
                proj += projection_inbreeding(partition, n_subsampling) * part_prob
            projection_matrix[allele_freq, :] = proj
        else:
            projection_matrix[allele_freq, :] = dadi.Numerics._cached_projection(n_subsampling, n_sequenced, allele_freq)
    
    return projection_matrix