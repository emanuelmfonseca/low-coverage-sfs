import numpy

def compute_cov_dist(data_dict, pop_ids):
    """
    Compute the coverage distribution for each population.

    Args:
    - data_dict (dict): A dictionary containing data entries.
    - pop_ids (list): A list of population identifiers for which coverage distribution is computed.

    Returns:
    - dict: A dictionary where keys are population identifiers, and values are arrays representing
            the coverage distribution for each population.
    
    Raises:
    - ValueError: If information about allelic depths for the reference and alternative alleles
                  is not found in the data dictionary.
    """
    try:
        # Dictionary comprehension to compute the coverage distribution
        coverage_distribution = {pop: numpy.array([
            *numpy.unique(numpy.concatenate([numpy.array(list(entry['coverage'][pop])) for entry in data_dict.values()]), return_counts=True)
        ]) for pop in pop_ids}
        
        # Normalize counts
        coverage_distribution = {pop: numpy.array([elements, counts / counts.sum()]) for pop, (elements, counts) in coverage_distribution.items()}
        
        return coverage_distribution
    except:
        raise ValueError("Information about allelic depths for the reference and alternative alleles not found in the data dictionary")