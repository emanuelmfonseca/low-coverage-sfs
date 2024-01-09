import numpy
from itertools import combinations

def projection_inbreeding(partition, k):
    """
    Calculate the distribution of inbreeding coefficients for a given partition.
    
    Args:
    - partition (iterable): A collection of indices representing the individuals in the partition.
    - k (int): The total number of individuals considered in each combination.
    
    Returns:
    - numpy.ndarray: An array representing the distribution of inbreeding coefficients. Each index in the array corresponds 
      to the total number of shared alleles in a combination, and the values represent the corresponding frequencies normalized 
      by the total number of combinations.
    ```
    """
    result = numpy.zeros_like(range(k+1))

    # Generate all possible combinations of partition indices
    partitions = list(combinations(partition, k//2))

    for p in partitions:
        # Calculate the sum for each partition and update the result array
        result[sum(p)] += 1

    return result/sum(result)