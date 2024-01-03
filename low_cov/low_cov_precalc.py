import numpy
from dadi.LowCoverage import probability_of_no_call_1D, probability_enough_individuals_covered, projection_matrix, calling_error_matrix, simulate_genotype_calling

def low_cov_precalc(nsub, nseq, cov_dist, sim_threshold=1e-2, Fx=0, nsim=1000):
    """
    Calculate transformation matrices for low-coverage calling model.
    
    nsub: Final sample size (in haplotypes)
    nseq: Sequenced sample size (in haplotypes)
    cov_dist: Coverage distribution (list of one array per population)
    sim_threshold: This method uses the probability an allele is not called
                   to switch between analytic and simulation-based methods.
                   Setting this threshold to 0 will always use simulations,
                   while setting it to 1 will always use analytics.
    nsim: For simulations, number of simulations per allele frequency combination
    """
    # As a lower bound on the probability that a allele with a given frequency is not called,
    # use the probability it is not called considering only the reads in each individual population.
    # We calculate them separately, then combine them into a single matrix
    prob_nocall_by_pop = [probability_of_no_call_1D(cov_dist_, nseq_, Fx_) for (cov_dist_, nseq_, Fx_) in zip(cov_dist.values(), nseq, Fx)]
    prob_nocall_ND = 1
    for apn_1D in prob_nocall_by_pop:
        prob_nocall_ND = numpy.multiply.outer(prob_nocall_ND, apn_1D)
    
    # Identify those entries that should be simulated, as opposed to analytically calculated
    use_sim_mat = prob_nocall_ND > sim_threshold
    
    ### For analytic calling model
    # Probability that enough individuals got at least one read to subsample
    # This doesn't depend on allele frequency
    prob_enough_covered = numpy.prod([probability_enough_individuals_covered(cov_dist_, nseq_, nsub_) for (cov_dist_, nseq_, nsub_) in zip(cov_dist.values(), nseq, nsub)])
    # Precalculate analytic projection and heterozygote error matrices
    proj_mats = [prob_enough_covered * projection_matrix(nseq_, nsub_, Fx_) for (nseq_, nsub_, Fx_) in zip(nseq, nsub, Fx)]
    heterr_mats = [calling_error_matrix(cov_dist_, nsub_, Fx_) for (cov_dist_, nsub_, Fx_) in zip(cov_dist.values(), nsub, Fx)]
    
    ### For simulation calling model
    # Indices of entries in afs where we should use simulation
    simulated_indices = numpy.argwhere(use_sim_mat)
    # Do the simulations (Note, arbitrarily using 1000 per entry here. Should think harder about tradeoffs.)
    sim_outputs = {tuple(af):simulate_genotype_calling(cov_dist, af, nseq, nsub, nsim, Fx) for af in simulated_indices}
    
    return prob_nocall_ND, use_sim_mat, proj_mats, heterr_mats, sim_outputs