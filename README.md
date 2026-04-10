# Behavioral Data-Driven Trajectory Generation for a Tower Crane

Supplementary MATLAB code and data for the paper:

> **Data-Driven Behavioral Trajectory Generation for a Tower Crane**  
> Iskandar Khemakhem, 2024

---

## Quick Start

```matlab
% From the MATLAB command window, navigate to this folder and run:
addpath(genpath(pwd))

% Then run any of the top-level scripts (see Workflow below).
```

**Dependencies:** MATLAB, CVX (http://cvxr.com/cvx), Signal Processing Toolbox.  
CVX must be installed and initialised (`cvx_setup`) before running scripts that call `simulate_dd_system` or `generate_trajectory`.  
The manuscript mentions MOSEK as a solver backend; the code is written at the CVX level and works with the default CVX solver as well.

---

## Repository Structure

```
.
├── data/
│   ├── crane_data/                  processed experimental crane trajectories
│   ├── sim_500/                     experimental validation set (20 sequences)
│   ├── Camera_Coordinates_June_Mean.mat  camera extrinsic calibration
│   ├── sim_grid_search.mat          precomputed best (threshold, lambda, p) — simulative study
│   └── trajopt_grid_search_long.mat precomputed best (mu, sigma, lambda) — trajectory generation
│
├── utils/                           core data-driven functions (add to MATLAB path)
│   ├── simulate_crane_one_step.m    one-step crane ODE integration (6-state model)
│   ├── simulate_model.m             full trajectory simulation from initial conditions
│   ├── generatePageMatrix.m         generate simulated data + block Hankel matrix
│   ├── pageMatrixFromData.m         load bundled crane data + build data matrix
│   ├── blkhank.m                    block Hankel matrix construction
│   ├── blkPage.m                    block Page matrix construction
│   ├── simulate_dd_system.m         nonparametric behavioral simulation (CVX)
│   ├── generate_trajectory.m        behavioral optimal trajectory generation (CVX)
│   ├── quantify_function.m          trajectory quality cost function
│   ├── organize_measured_trajectory.m  preprocess one experimental result file
│   ├── pixels_to_world.m            pixel-to-angle conversion using camera calibration
│   ├── linear_kalman.m              Kalman filter for boom angle/velocity estimation
│   ├── calculateAcceleration.m      numerical differentiation of velocity profiles
│   ├── generate_random_input.m      smooth random input generation
│   ├── custom_spaced_values.m       skewed parameter grid generation
│   ├── time_to_reach_target.m       settling-time computation
│   └── vec.m                        column-vectorisation helper
│
├── model_based/                     model-based benchmark (comparison only)
│   ├── MAIN.m                       entry point for model-based trajectory optimisation
│   ├── Verification.m               crane load-sway ODE (4-state, used by MAIN.m)
│   ├── Verification_Model.m         simplified crane ODE variant
│   ├── nonlcon.m                    nonlinear constraints for fmincon
│   ├── obj_function.m               objective function for fmincon
│   ├── Create_Data_Mat_Files.m      export optimised trajectory to .mat
│   └── Compare_Results.m            plot data-driven vs. model-based results
│
├── simulation_validation.m          validate nonparametric simulation on simulated data
├── grid_search_nonparametricSim_simulative.m  tune (threshold, lambda, p) — simulated data
├── real_data_grid_search_simulation.m         tune (threshold, lambda, p) — experimental data
├── traj_generation_grid_search.m              tune (mu, sigma, lambda) — experimental data
└── real_data_grid_search_opt_traj.m           generate optimal trajectory — experimental data
```

---

## Workflow

Scripts are designed to be run from the repository root with `addpath(genpath(pwd))` active.

### Simulative pipeline (no real crane required)

| Step | Script | Output |
|------|--------|--------|
| 1 | `simulation_validation.m` | Validates nonparametric simulation against the model |
| 2 | `grid_search_nonparametricSim_simulative.m` | Tunes `(threshold, lambda, p)`; saves `data/sim_grid_search.mat` |

### Experimental pipeline (uses bundled crane data)

| Step | Script | Output |
|------|--------|--------|
| 1 | `real_data_grid_search_simulation.m` | Tunes `(threshold, lambda, p)` on real data; saves `data/sim_grid_search_real_data_500.mat` |
| 2 | `traj_generation_grid_search.m` | Tunes `(mu, sigma, lambda)`; saves `data/trajopt_grid_search.mat` |
| 3 | `real_data_grid_search_opt_traj.m` | Generates the optimal trajectory using precomputed best parameters |

Steps 2 and 3 of the experimental pipeline can be run immediately with the precomputed parameter files already included in `data/`.

---

## Data

### `data/crane_data/`
29 processed experimental trajectories from the BehavioralCrane laboratory setup.  
Each sequence consists of:
- `input_seq_NN.mat` — variable `acceleration` (boom angular acceleration `ddtheta4`, `[rad/s²]`)
- `results_seq_NN.mat` — variable `yd_i` (`[theta1, theta2, theta4, dtheta4]`, `[rad, rad, rad, rad/s]`)

Loaded by `utils/pageMatrixFromData.m`.

### `data/sim_500/`
20 experimental validation sequences used by `real_data_grid_search_simulation.m`:
- `input_seq_sim_500_N.mat` — variables `acceleration` and `velocity`
- `results_seq_sim_500_N.mat` — variables `raw_vision_data` and `vision_data`

### `data/Camera_Coordinates_June_Mean.mat`
Extrinsic calibration of the RGB camera used for load-position measurement.  
Contains `A_WP` (rotation matrix, world←camera) and `W_r_OP` (camera position in world frame),
obtained by averaging three PnP solutions at fixed camera mounting.  
Loaded automatically by `utils/pixels_to_world.m`.

### `data/sim_grid_search.mat`
Precomputed best `(threshold, lambda, p)` from the simulative nonparametric-simulation study.

### `data/trajopt_grid_search_long.mat`
Precomputed best `(mu, sigma, lambda)` for trajectory generation (long grid search).

---

## Camera Setup

Load-position measurements use a single standard RGB camera (1280×720, 20 fps) mounted
above the crane. The boom carries a coloured marker; pixel coordinates are extracted by
colour-blob detection and converted to crane angles `(theta1, theta2)` via
`utils/pixels_to_world.m` using the bundled PnP calibration.  
Boom angle `theta4` and velocity `dtheta4` are recovered from encoder readings and
smoothed with the linear Kalman filter in `utils/linear_kalman.m`.

---

## Note on `q = 5`

The stacked trajectory sample is

```
w(i) = [ddot(theta_4), theta_1, theta_2, theta_4, dot(theta_4)]^T
```

so `q = 5`. The first channel (`ddtheta4`) is included as an input because the behavioral
framework treats inputs and outputs symmetrically. `theta4` and `dtheta4` are both retained
so that boundary conditions on the boom position and velocity can be imposed directly in the
trajectory-generation problem.

---

## Note on `model_based/`

The `model_based/` folder contains the model-based benchmark used for comparison in the paper.
`Verification.m` implements the 4-state load-sway ODE driven by prescribed boom velocity and
acceleration; it is called from `MAIN.m` via `fmincon`. This is independent of the data-driven
pipeline: the data-driven scripts use `utils/simulate_crane_one_step.m`, which embeds the full
6-state crane dynamics (including integrated boom motion) with hardcoded physical parameters.
