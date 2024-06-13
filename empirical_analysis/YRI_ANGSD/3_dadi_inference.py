#!/usr/bin/env python

#SBATCH --job-name=dadi_inference_YRI_ANGSD
#SBATCH --output=outfiles/%x-%j.out
#SBATCH --error=outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=02:00:00

# micromamba activate low-coverage when submit job
# conda activate low_coverage has an error because of old python/numpy version

import os
import dadi
import nlopt
import numpy as np
from dadi.LowCoverage.LowCoverage import *

fs_dir="angsd_results"
result_dir="inference_two_epoch"

# make results dir
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

fn_list=[
        '3x.sfs.mle',
        '30x.sfs.mle',
        '5x.sfs.mle',
        '10x.sfs.mle',
        ]
        
for fn in fn_list:
    # load ANGSD fs
    with open(os.path.join(fs_dir, fn+'.txt')) as fh:
        fs_txt = fh.readline().strip()
        fs = [float(entry) for entry in fs_txt.split()]
        data_fs = dadi.Spectrum(fs)
    
    # Retrive the sample sizes from the data
    ns = data_fs.sample_sizes

    # Define the grid points based on the sample size.
    # For smaller data (largest sample size is about <100) [ns+20, ns+30, ns+40] is a good starting point.
    # For larger data (largest sample size is about >=100) or for heavily down projected data [ns+100, ns+110, ns+120] is a good starting point.
    pts_l = [max(ns)+120, max(ns)+130, max(ns)+140]
    
    # Define starting parameters
    params = [1, 0.01, 0.01]

    # Define boundaries of optimization.
    # It is a good idea to have boundaries to avoid optimization
    # from trying parameter sets that are time consuming without
    # nessicarily being correct.
    # If optimization infers parameters very close to the boundaries, we should increase them.
    lower_bounds = [1e-2, 1e-3, 1e-3]
    upper_bounds = [1e2, 3, 1]
    
    # init demographic model, regular (not applying low coverage)
    demo_model = dadi.Numerics.make_anc_state_misid_func(dadi.Demographics1D.two_epoch)
    # demo_model = dadi.Numerics.make_anc_state_misid_func(dadi.Demographics1D.growth)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    n_rep=100
    # first optimize normally (no low_cov correction)
    for _ in range(n_rep):
        # Create or append to an file to store optimization results
        try:
          fid = open(os.path.join(result_dir, fn+'_chr20_1d_demo_fits.txt'),'a')
        except:
          fid = open(os.path.join(result_dir, fn+'_chr20_1d_demo_fits.txt'),'w')
          
        p0 = dadi.Misc.perturb_params(params, 
                                     upper_bound=upper_bounds,
                                     lower_bound=lower_bounds)
    
        popt, ll_model = dadi.Inference.opt(p0, data_fs, 
                                     demo_model_ex, pts_l,
                                     lower_bound=lower_bounds,
                                     upper_bound=upper_bounds,
                                     algorithm=nlopt.LN_BOBYQA,
                                     maxeval=600, verbose=0)
    
        # Calculate the theta
        model_fs = demo_model_ex(popt, ns, pts_l)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, data_fs)
    
        # Write results to fid
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res])+'\n')
    
    fid.close()

    # sort the output
    infile_path=os.path.join(result_dir, fn+'_chr20_1d_demo_fits.txt')
    outfile_path=os.path.join(result_dir, 'sorted_' + fn + '_chr20_1d_demo_fits.txt')
    myCmd = f'sort {infile_path} > {outfile_path}'
    os.system(myCmd)

