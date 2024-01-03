import dadi
import numpy
import math

def part_inbreeding_probability(parts, Fx):
    """
    Calculate genotype partition probabilities under inbreeding.

    Parameters:
    - parts (list of lists): List of genotype partitions.
    - Fx (float): Inbreeding coefficient.

    Returns:
    - numpy.array: Normalized genotype partition probabilities.
    """
    part_prob = numpy.array([])
    for part in parts:
        if sum(part) != 0 and sum(part) != 2*len(part):
            p = (2 * part.count(2) + part.count(1))/(2*len(part))
            alpha = p*((1.0-Fx)/Fx)
            beta  = (1.0-p)*((1.0-Fx)/Fx)
            
            p00, p01, p11 = numpy.exp([dadi.Numerics.BetaBinomln(_,2,alpha,beta) for _ in range(2+1)])
            n, n00, n01, n11 = len(part), part.count(0), part.count(1), part.count(2)
            
            part_prob = numpy.append(part_prob, (math.factorial(n) / (math.factorial(n00) * math.factorial(n01) * math.factorial(n11))) * (p00 ** n00) * (p01 ** n01) * (p11 ** n11))
        else:
            part_prob = numpy.append(part_prob, 1)
        
    return part_prob / sum(part_prob)