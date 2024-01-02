# Optimizing demographic inference with dadi using low-coverage sequencing data

## Overview

This repository contains an implementation of a method for correcting low-coverage distortion in demographic inference. The method is designed to address challenges arising from insufficient data coverage in demographic studies, enhancing the accuracy and reliability of the results. Within this repository, you will find code, data, and results encompassing all simulations and empirical data analyses presented in the manuscript 'Inferring beyond the blur: Unveiling population histories with corrected low-coverage genomes.

## Preprint

Fonseca, E.M., Tran, L., Mendonza, H., Gutenkunst, R.N. Beyond the blur: unveiling population histories with corrected low-coverage genomes (coming soon).

## What does this repository contain?

 - The `/low-cov` directory contains the entire codebase used to implementing the method designed to correct bias arising from low-coverage sequencing.
 
 -  The `/simulations` directory contains the code for generating and running simulations at various coverage levels. Moreover, it includes the resulting site frequency spectra from these simulations

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

