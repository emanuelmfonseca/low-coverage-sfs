# SFS-based model to optimize demographic inference by addressing low-coverage sequencing bias

## Overview

This repository contains an implementation of a method for correcting low-coverage distortion in demographic inference. The method is designed to address challenges arising from insufficient data coverage in demographic studies, enhancing the accuracy and reliability of the results. Within this repository, you will find code, data, and results encompassing all simulations and empirical data analyses presented in the manuscript 'Inferring beyond the blur: Unveiling population histories with corrected low-coverage genomes.

## Prepring

Fonseca, E.M., Tran, L., Mendonza, H., Gutenkunst, R.N. Inferring beyond the blur: Unveiling population histories with corrected low-coverage genomes (coming soon).

## Getting Started

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
make_low_cov_func(demo_model, data_dict, pop_ids, nseq, nsub, sim_threshold, inbreeding), where:

- demo_model: specified demographic model in dadi
- data_dict: a data dictionary comprising information extracted from a VCF file
- pop_ids: population names to be analyzed
- nseq: total number of samples for a given population
- nsub: subsampled number of samples for a given population
- sim_threshold: This method switches between analytic and simulation-based methods. Setting this threshold to 0 will always use simulations, while setting it to 1 will always use analytics. Values in between indicate that simulations will be employed for thresholds below that value.
- inbreeding: should inbreeding be included? [True/False]

```
