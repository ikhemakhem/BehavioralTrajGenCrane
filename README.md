<<<<<<< HEAD
# Complementary Material

This folder contains the supplementary MATLAB code and data for the behavioral data-driven trajectory-generation method.


## Folder Structure

- `data/`
  Bundled data used by the scripts.
- `utils/`
  MATLAB helper functions for data processing, behavioral simulation, and trajectory generation.
- `utils/` 
   scrit used for the comparison with model-based benchmark.
- `grid_search_nonparametricSim_simulative.m`
  Simulative nonparametric-simulation study and tuning of `(delta, lambda, nu)`.
- `traj_generation_grid_search.m`
  Simulative behavioral trajectory generation and tuning of `(mu, sigma, lambda)`.
- `real_data_grid_search_simulation.m`
  Experimental nonparametric model-validation workflow.
- `real_data_grid_search_opt_traj.m`
  Experimental behavioral trajectory-generation workflow.

## What Is Included

This release contains the core behavioral / nonparametric method used in the paper:

- data-matrix construction via block Hankel / Page matrices,
- sparse missing-data recovery / nonparametric simulation,
- QR-based column selection and SVD truncation,
- grid-search tuning of simulation and trajectory-generation hyperparameters,
- direct behavioral optimal trajectory generation with smoothness and boundary constraints,
- experimental preprocessing utilities such as pixel-to-world conversion and Kalman filtering.


## Data

The bundled data are organized as follows:

- `data/crane_data/`
  Processed experimental trajectories used by `utils/pageMatrixFromData.m`.
- `data/sim_500/`
  Simulative validation trajectories used for the experimental nonparametric-simulation study.
- `data/sim_grid_search.mat`
  Precomputed best parameters for the simulative nonparametric-simulation study.
- `data/trajopt_grid_search_long.mat`
  Precomputed best parameters for trajectory generation.

`utils/pageMatrixFromData.m` is the portable loader for the bundled processed crane data.

`utils/pageMatrixData.m` is retained for the original raw-data workflow and still assumes a machine-specific local raw-data directory.

`real_data_grid_search_simulation.m` uses the `data/sim_500/` validation set. If you run this script directly, make sure those files are on the MATLAB path or adapt the load paths accordingly.

## Note on `q = 5`

Some scripts use the stacked trajectory sample

`w(i) = [ddot(theta_4), theta_1, theta_2, theta_4, dot(theta_4)]^T`

so `q = 5`.

This is intentional. The extra channel is the boom angular acceleration. It is included to keep the behavioral formulation compatible with the simulative setting while still retaining the boom velocity in the trajectory vector.

So `q = 5` should be understood as an implementation detail of the stacked trajectory representation, not as a change of method.

## Dependencies

The code relies on:

- MATLAB,
- CVX,
- Signal Processing Toolbox.

The manuscript mentions MOSEK. In practice, MOSEK can be used as the CVX solver backend, but the code is written at the CVX level.

