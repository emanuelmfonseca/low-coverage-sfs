# micromamba activate low-coverage
import os
import glob
import dadi
import numpy as np
from matplotlib import pyplot as plt

if not os.path.exists('inference_two_epoch/plots'):
    os.makedirs('inference_two_epoch/plots')

depth_list=['3x',
            '5x',
            '10x',
            '30x',
            ]

for depth in depth_list:
    # data fs
    data_fs = dadi.Spectrum.from_file(f'inference_two_epoch/{depth}_fs_subsampled_32')
    
    # all model fs
    model_fs_list = sorted(glob.glob(f"inference_two_epoch/{depth}_fs_subsample*model_[0-9]"))
    low_cov_model_fs_list = sorted(glob.glob(f"inference_two_epoch/{depth}_fs_subsample*low_cov_[0-9]"))
    
    # pick the best fs
    model_fs = dadi.Spectrum.from_file(model_fs_list[0])
    low_cov_model_fs = dadi.Spectrum.from_file(low_cov_model_fs_list[0])
    
    # plot for no correction vs. low coverage correction
    fig = plt.figure()
    dadi.Plotting.plot_1d_comp_multinom(model_fs, data_fs)
    plt.title(f'{depth}')
    plt.savefig(f'inference_two_epoch/plots/{depth}.png', dpi=150)
    
    fig = plt.figure()
    dadi.Plotting.plot_1d_comp_multinom(low_cov_model_fs, data_fs)
    plt.title(f'{depth}_low_cov')
    plt.savefig(f'inference_two_epoch/plots/{depth}_low_cov.png', dpi=150)
