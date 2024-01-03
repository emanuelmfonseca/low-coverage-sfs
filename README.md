# Optimizing demographic inference with dadi using low-coverage sequencing data

## Overview

This repository contains an implementation of a method for correcting low-coverage distortion in demographic inference. The method is designed to address challenges arising from insufficient data coverage in demographic studies, enhancing the accuracy and reliability of the results. Within this repository, you will find code, data, and results encompassing all simulations and empirical data analyses presented in the manuscript 'Inferring beyond the blur: Unveiling population histories with corrected low-coverage genomes.

## Preprint

Fonseca, E.M., Tran, L., Mendonza, H., Gutenkunst, R.N. Beyond the blur: unveiling population histories with corrected low-coverage genomes (coming soon).

## What does this repository contain?

 - The `/low-cov` directory contains the entire codebase used to implementing the method designed to correct bias arising from low-coverage sequencing.
 
 -  The `/simulations` directory contains the code for generating and running simulations at various coverage levels. Moreover, it includes the resulting site frequency spectra from these simulations.

 -  The `/demography` directory contains the code for generating and running demographic inference under various coverage levels.

### Prerequisites

- Python 3
- dadi (https://dadi.readthedocs.io/en/latest/)
- numpy
- math
- scipy
- itertools
- warnings

### Usage

```python
make_low_cov_func(demo_model, data_dict, pop_ids, nseq, nsub, sim_threshold, inbreeding)
```

- demo_model: specified demographic model in dadi
- data_dict: a data dictionary comprising information extracted from a VCF file
- pop_ids: population names to be analyzed
- nseq: total number of samples for a given population
- nsub: subsampled number of samples for a given population
- sim_threshold: This method switches between analytic and simulation-based methods. Setting this threshold to 0 will always use simulations, while setting it to 1 will always use analytics. Values in between indicate that simulations will be employed for thresholds below that value.
- inbreeding: should inbreeding be included? [True/False]

### Example

#### Import necessary libraries
```python
import dadi
from dadi.LowCoverage import LowCoverage
import nlopt
```

#### Paths to input data files
```python
datafile = '/simulations/simulated_datasets/1D_exp_growth_two_epochs_nu1_10_T1_0.1/heterogeneous_coverage/coverage_3/VCF_files_gatk/gatk_Replicate_1_filtered.vcf'
popfile = '/simulations/simulated_datasets/1D_exp_growth_two_epochs_nu1_10_T1_0.1/heterogeneous_coverage/coverage_3/VCF_files_gatk/popfile.txt'
```

#### Define data parameters
```python
pop_id = 'pop1'  # Population name
nseq = 40  # Total number of sequenced chromosomes
nsub = 32  # Number of chromosomes to be subsampled
ss = {pop_id: nsub // 2}  # Subsampling dictionary
```

#### Create data dictionary from the VCF file
```python
data_dict = dadi.Misc.make_data_dict_vcf(datafile, popfile, subsample=ss)
```

#### Assign outgroup information to each genomic position
```python
for chrom_pos in data_dict:
    data_dict[chrom_pos]['outgroup_allele'] = data_dict[chrom_pos]['segregating'][0]
    data_dict[chrom_pos]['outgroup_context'] = data_dict[chrom_pos]['segregating'][0]
```

#### Generate the Site Frequency Spectrum (SFS)
```python
data_fs = dadi.Spectrum.from_data_dict(data_dict, [pop_id], [nsub])
```

#### Define the demographic model (exponential growth)
```python
demo_model_ex = dadi.Numerics.make_extrap_func(dadi.Demographics1D.growth)
```

#### Wrap the demographic model with low coverage model
```python
demo_model_ex = LowCoverage.make_low_cov_func(demo_model_ex, data_dict, data_fs.pop_ids, [nseq], [nsub], sim_threshold=1e-2, inbreeding=False)
```

#### Define demographic parameters
```python
param_names = ['nu1', 'T']  # Growth rate and time since growth started
lower_bounds = [1e-05, 1e-05]  # Lower bounds for parameters
upper_bounds = [50, 1]  # Upper bounds for parameters
params = [1, 0.05]  # Initial parameter values
grids = [50, 60, 70]  # Grid points for optimization
```

#### Perturb the initial parameters
```python
p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)
```

#### Run optimization for parameter inference
```python
popt, ll_model = dadi.Inference.opt(p0, data_fs, demo_model_ex, grids,
                                    lower_bound=lower_bounds,
                                    upper_bound=upper_bounds,
                                    algorithm=nlopt.LN_COBYLA,
                                    maxeval=1000, verbose=0)
```

#### Print the results
```python
print('True parameters are: nu1 = 10 and T = 0.1')
print(f'The inferred parameters were: nu1 = {popt[0]:.2f} and T = {popt[1]:.2f}')
```

### Questions

For any questions, please contact Emanuel M. Fonseca (emanuelmfonseca@arizona.email.edu or emanuelmfonseca@gmail.com) or Ryan Gutenkunst (rgutenk@arizona.edu).

