import os
import sys
import dadi
import nlopt
from matplotlib import pyplot as plt
from dadi.LowCoverage.LowCoverage import *

vcf_dir="gatk_called_vcf"
result_dir="inference_split_mig"

pop_id = ['YRI', 'CEU']  # Population name
nseq = [20, 20]  # Total number of sequenced chromosomes
nsub = [16, 16]  # Number of chromosomes to be subsampled
ss = {'YRI': int(nsub[0]//2), 'CEU': int(nsub[1]//2)} # Subsampling dictionary

prefix='aa_annotated_CEU_YRI'
suffix='GRCh38_chr20_gatk_filtered.vcf.gz'

depth = sys.argv[1]

# make subsampled data dict from vcf
dd = dadi.Misc.make_data_dict_vcf(os.path.join(vcf_dir, f'{prefix}_{depth}_{suffix}'), 'popfile.txt', subsample=ss)

# make dadi fs from vcf
data_fs = dadi.Spectrum.from_data_dict(dd, pop_id, nsub)

# Save fs to disk
data_fs.to_file(os.path.join(result_dir, f'{depth}_fs_subsampled_16'))

# # plot fs
# try:
#     dadi.Plotting.plot_single_2d_sfs(data_fs)
#     plt.savefig(os.path.join(result_dir, f'{depth}_fs_subsampled_16.png'), dpi=150)
# except ValueError:
#     print("Invalid vmin or vmax for plotting")

# Retrive the sample sizes from the data
ns = data_fs.sample_sizes

# Define the grid points based on the sample size.
# For smaller data (largest sample size is about <100) [ns+20, ns+30, ns+40] is a good starting point.
# For larger data (largest sample size is about >=100) or for heavily down projected data [ns+100, ns+110, ns+120] is a good starting point.
pts_l = [max(ns)+120, max(ns)+130, max(ns)+140]

# Define starting parameters
params = [1, 1, 0.05, 1, 0.01]

# Define boundaries of optimization.
# It is a good idea to have boundaries to avoid optimization
# from trying parameter sets that are time consuming without
# nessicarily being correct.
# If optimization infers parameters very close to the boundaries, we should increase them.
lower_bounds = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3]
upper_bounds = [1e2, 1e2, 3, 10, 1]

# init demographic model, regular (not applying low coverage)
demo_model = dadi.Numerics.make_anc_state_misid_func(dadi.Demographics2D.split_mig)
demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

n_rep=100
# first optimize normally (no low_cov correction)
for _ in range(n_rep):
    # Create or append to an file to store optimization results
    try:
      fid = open(os.path.join(result_dir, f'{depth}_2d_demo_fits.txt'),'a')
    except:
      fid = open(os.path.join(result_dir, f'{depth}_2d_demo_fits.txt'),'w')
      
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
infile_path=os.path.join(result_dir, f'{depth}_2d_demo_fits.txt')
outfile_path=os.path.join(result_dir, f'sorted_{depth}_2d_demo_fits.txt')
myCmd = f'sort {infile_path} > {outfile_path}'
os.system(myCmd)

# save the top three model FS for plotting later
# readin the sorted inference file
print(f'Top optimizations for {outfile_path}:')
with open(outfile_path) as fh:
    head = [next(fh).strip().split() for _ in range(3)]
    for i, inference in enumerate(head):
        inference_num = [round(float(num),3) for num in inference]
        print(*inference_num, sep='\t')
        
        # get popt for reconstructing the fs
        popt = inference_num[1:-1] # ignore ll and theta
        model_fs = demo_model_ex(popt, ns, pts_l)
        model_fs.to_file(os.path.join(result_dir, f'{depth}_fs_subsampled_16_model_{i}'))


# low coverage version
demo_model_ex = make_low_cov_func(demo_model_ex, dd, data_fs.pop_ids, nseq, nsub)

# optimize with low_coverage
for _ in range(n_rep):
    # Create or append to an file to store optimization results
    try:
      fid = open(os.path.join(result_dir, f'{depth}_2d_demo_fits_low_cov.txt'),'a')
    except:
      fid = open(os.path.join(result_dir, f'{depth}_2d_demo_fits_low_cov.txt'),'w')
      
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
infile_path=os.path.join(result_dir, f'{depth}_2d_demo_fits_low_cov.txt')
outfile_path=os.path.join(result_dir, f'sorted_{depth}_2d_demo_fits_low_cov.txt')
myCmd = f'sort {infile_path} > {outfile_path}'
os.system(myCmd)

# save the top three model FS for plotting later
# readin the sorted inference file
print(f'Top optimizations for {outfile_path}:')
with open(outfile_path) as fh:
    head = [next(fh).strip().split() for _ in range(3)]
    for i, inference in enumerate(head):
        inference_num = [round(float(num),3) for num in inference]
        print(*inference_num, sep='\t')
        
        # get popt for reconstructing the fs
        popt = inference_num[1:-1] # ignore ll and theta
        model_fs = demo_model_ex(popt, ns, pts_l)
        model_fs.to_file(os.path.join(result_dir, f'{depth}_fs_subsampled_16_model_low_cov_{i}'))
        