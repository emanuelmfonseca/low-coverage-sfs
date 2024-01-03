import numpy
from itertools import combinations


def projection_inbreeding(partition, k):
    result = numpy.zeros_like(range(k+1))

    # Generate all possible combinations of partition indices
    partitions = list(combinations(partition, k//2))

    for p in partitions:
        # Calculate the sum for each partition and update the result array
        result[sum(p)] += 1

    return result/sum(result)