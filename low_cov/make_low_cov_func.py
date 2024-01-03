import numpy
from dadi.LowCoverage import compute_cov_dist, low_cov_precalc 

def make_low_cov_func(func, dd, pop_ids, nseq, nsub, sim_threshold=1e-2, inbreeding=False):
    """
    Generage a version of func accounting for low coverage distortion
    
    func: The frequency spectrum generating function to which distortion will
          be applied. It is assumed that the second argument is the final
          sample size (in haplotypes).
    vcf_file: The original data VCF file, from which coverage information will
              extracted.
    popfile: Name of file containing the population assignments for
             each sample in the VCF. If a sample in the VCF file does
             not have a corresponding entry in this file, it will be
             skipped. See dadi.Misc._get_popinfo for information on how this
             file must be formatted.
    nseq: Number of sequenced individuals. # XXX: We should just extract this when we parse the VCF file.
    sim_threshold: This method uses the probability an allele is not called
                   to switch between analytic and simulation-based methods.
                   Setting this threshold to 0 will always use simulations,
                   while setting it to 1 will always use analytics.
    """
    # XXX: Should automatically extract nseq from the vcf file
    cov_dist = compute_cov_dist(dd, pop_ids)
    
    # Used to cache matrices used for low-coverage transformation
    precalc_cache = {}
    
    def lowcov_func(*args, **kwargs):
        ns = args[1]
        Fx = args[0][-len(nseq):] if inbreeding else [0] * len(ns)
        new_args = [args[0]] + [nseq] + list(args[2:])
        model = func(*new_args, **kwargs)
        if model.folded:
            raise ValueError('Low coverage model not tested for folded model spectra yet.')
        
        if tuple(nsub) not in precalc_cache:
            precalc_cache[tuple(nsub)] = low_cov_precalc(nsub, nseq, cov_dist, sim_threshold, Fx)
        prob_nocall_ND, use_sim_mat, proj_mats, heterr_mats, sim_outputs = precalc_cache[tuple(nsub)]
        # First, transform entries we do analytically. We zero out the entries
        # we'll simulate, since we'll handle their contribution later.
        analytic = model * (1-use_sim_mat)
        # Account for sites that aren't called 
        analytic *= (1-prob_nocall_ND)
        # Apply projection and heterozygote error transformations
        for pop_ii, (proj_mat, heterr_mat) in enumerate(zip(proj_mats, heterr_mats)):
            analytic = analytic.swapaxes(pop_ii, -1) # Swap axes to use efficient matrix multiplication
            analytic = analytic.dot(proj_mat)
            analytic = analytic.dot(heterr_mat)
            analytic = analytic.swapaxes(pop_ii, -1)
        
        # Use the simulated outputs
        simulated = numpy.sum([model[af]*output for (af, output) in sim_outputs.items()], axis=0)
        
        output = analytic + simulated
        # Not sure why the folding status got undefined in the above manipulations...
        output.folded = model.folded
        output.extrap_x = model.extrap_x
        
        return output
    lowcov_func.__name__ = func.__name__ + '_lowcov'
    lowcov_func.__doc__ = func.__doc__
    return lowcov_func