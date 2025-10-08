# MIMO_HMC_detector

Simulation codes for MIMO signal detection using Hamiltonian Monte Carlo (HMC).

Paper: [Hamiltonian Monte Carlo-Based Near-Optimal MIMO Signal Detection](https://arxiv.org/abs/2412.02391)

## Overview

This repository implements near-optimal MIMO signal detection by transforming the discrete detection problem into a continuous framework and leveraging the HMC algorithm via Stan.

## Requirements

- R >= 4.0
- rstan
- ggplot2

## Installation

```r
install.packages("rstan")
install.packages("ggplot2")
```

## Usage

```r
source("main.R")
# Run simulation
```

## Citation

```bibtex
@article{hagiwara2024hamiltonian,
  title={Hamiltonian Monte Carlo-Based Near-Optimal MIMO Signal Detection},
  author={Hagiwara, Junichiro et al.},
  journal={arXiv preprint arXiv:2412.02391},
  year={2024}
}
```

## License

MIT License

## Status

Under development for IEEE Transactions submission.
