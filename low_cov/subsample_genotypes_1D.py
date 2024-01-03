import numpy

rng = numpy.random.default_rng()

def subsample_genotypes_1D(genotype_calls, n_subsampling):
    """
    Subsample genotypes to create a smaller dataset with a specified number of haplotypes.
    
    Args:
        genotype_calls (numpy.ndarray): Genotype array with 99 assumed to represent missing data.
        n_subsampling (int): Number of haplotypes in the final subsample.
    
    Returns:
        numpy.ndarray: Subsampled genotype data.
    """
    # Handle a special case where the input genotype_calls is empty
    if len(genotype_calls) == 0:
        return genotype_calls[:, :n_subsampling // 2]
    
    # Count the number of called individuals for each locus
    n_called = numpy.count_nonzero(genotype_calls != 99, axis=1)
    
    # Sort the genotypes within each locus, placing uncalled genotypes at the end
    sorted_genotype_calls = numpy.sort(genotype_calls, axis=1)
    
    subsampled_data = []
    # Iterate through unique counts of called individuals
    for calls in numpy.sort(numpy.unique(n_called)):
        if calls < n_subsampling // 2:
            continue  # Skip loci without enough calls
        
        # Extract loci with exactly 'calls' called individuals and keep only the called genotypes
        loci_with_calls = sorted_genotype_calls[n_called == calls][:, :calls]
        
        # Permute the order of genotypes within each locus
        permuted_loci = rng.permuted(loci_with_calls, axis=1)
        
        # Take the first 'n_subsampling // 2' genotypes from each (permuted) locus
        subsampled_data.append(permuted_loci[:, :n_subsampling // 2])
    
    # Concatenate the subsamples from each group of called individuals
    return numpy.concatenate(subsampled_data)
