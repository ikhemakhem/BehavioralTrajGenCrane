# Behavioral Data-Driven Trajectory Generation for a Tower Crane

Supplementary MATLAB code and data for the paper:

> **Data-Driven Behavioral Trajectory Generation for a Tower Crane**  
> Iskandar Khemakhem, 2024

---

## Quick Start

```matlab
% From the MATLAB command window, navigate to the repository root and run:
addpath(genpath(pwd))

% Then run any top-level script (see Workflow below).
```

**Dependencies:** MATLAB, CVX (http://cvxr.com/cvx), Signal Processing Toolbox.  
CVX must be installed and initialised (`cvx_setup`) before running any script
that calls `simulate_dd_system` or `generate_trajectory`.  
The manuscript uses MOSEK as a solver backend; the code is written at the CVX
level and also works with the default CVX solver.

---

## Repository Structure

```
.
├── data/
│   ├── training_data/               29 experimental sequences — used to build the data matrix
│   ├── validation_data/             20 experimental sequences — used to evaluate predictions
│   ├── camera_coordinates.mat       camera PnP calibration (A_WP, W_r_OP)
│   ├── opt_hyperparams_nonparamSim.mat  precomputed best (threshold, lambda, p)
│   └── opt_hyperparams_trajOpt.mat      precomputed best (mu, sigma, lambda)
│
├── utils/                           core data-driven functions (add to MATLAB path)
│   ├── simulate_crane_one_step.m    one-step crane ODE integration (6-state model)
│   ├── simulate_model.m             full trajectory simulation from initial conditions
│   ├── generatePageMatrix.m         generate model data + block Hankel matrix
│   ├── pageMatrixFromData.m         load training data + build data matrix
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
├── simulation_validation.m          quick sanity check: nonparametric vs. model
│
├── -- Grid search (standalone, long-running) --
├── grid_search_nonparametric_sim.m  tune (threshold, lambda, p) — model-generated data
├── grid_search_nonparametric_real.m tune (threshold, lambda, p) — experimental data
├── grid_search_trajectory.m         tune (mu, sigma, lambda) — trajectory optimisation
│
└── -- Execution scripts (use precomputed or manual parameters) --
    ├── run_nonparametric_simulation.m   validate simulation on all 20 validation sequences
    └── run_trajectory_generation.m      generate optimal trajectory + results figure
```

---

## Data

Both data folders contain **real experimental recordings** from the physical
BehavioralCrane laboratory setup.  They differ in their role:

| Folder | Sequences | Role |
|--------|-----------|------|
| `data/training_data/` | 29 | Build the Hankel/Page data matrix H_w (the system representation used by all data-driven methods) |
| `data/validation_data/` | 20 | Evaluate predictions — sequences 1–15 serve as the grid-search objective; all 20 are used in the final validation plots |

### `data/training_data/`
Each sequence consists of:
- `input_seq_NN.mat` — variable `acceleration` (boom angular acceleration `ddtheta4`, `[rad/s²]`)
- `results_seq_NN.mat` — variable `yd_i` (`[theta1, theta2, theta4, dtheta4]`, shape `T×4`, `[rad]`)

Loaded automatically by `utils/pageMatrixFromData.m`.

### `data/validation_data/`
Each sequence consists of:
- `input_seq_sim_500_N.mat` — variables `acceleration` and `velocity`
- `results_seq_sim_500_N.mat` — variables `raw_vision_data` and `vision_data`

### `data/camera_coordinates.mat`
Extrinsic calibration of the single standard RGB camera (1280×720, 20 fps)
mounted above the crane.  Contains `A_WP` (rotation matrix, world←camera)
and `W_r_OP` (camera position in world frame), obtained by averaging three
PnP solutions at a fixed camera mounting.  
Loaded automatically by `utils/pixels_to_world.m`.

### `data/opt_hyperparams_nonparamSim.mat`
Precomputed best `(threshold, lambda, p)` for the nonparametric simulation,
obtained from `grid_search_nonparametric_sim.m`.

### `data/opt_hyperparams_trajOpt.mat`
Precomputed best `(mu, sigma, lambda)` for trajectory generation,
obtained from an extensive run of `grid_search_trajectory.m`.

---

## Workflow

Run scripts from the repository root with `addpath(genpath(pwd))` active.

### Simulative pipeline (no real crane required)

| Step | Script | Output |
|------|--------|--------|
| 1 | `simulation_validation.m` | Validates the nonparametric framework against the ODE model |
| 2 | `grid_search_nonparametric_sim.m` | Tunes `(threshold, lambda, p)`; saves `data/opt_hyperparams_nonparamSim.mat` |

### Experimental pipeline (uses bundled crane data)

Each execution script runs immediately with the precomputed parameter files
already included in `data/`.  Set `USE_PRECOMPUTED_* = false` at the top of
the script to supply your own hyperparameter values instead.

| Step | Script | Output |
|------|--------|--------|
| 1 | `grid_search_nonparametric_real.m` | Tunes `(threshold, lambda, p)` on real data; saves `data/opt_hyperparams_nonparamSim_real.mat` |
| 2 | `run_nonparametric_simulation.m` | Validates simulation on all 20 validation sequences (plots prediction vs. measurement) |
| 3 | `grid_search_trajectory.m` | Tunes `(mu, sigma, lambda)`; saves `data/opt_hyperparams_trajOpt.mat` |
| 4 | `run_trajectory_generation.m` | Generates optimal trajectory + results figure; saves `data/opt_traj_result.mat` |

Steps 2 and 4 can be run immediately using the precomputed files in `data/`.

---

## Camera Setup

Load-position measurements use a single standard RGB camera (1280×720, 20 fps)
mounted above the crane.  The boom carries a coloured marker; pixel coordinates
are extracted by colour-blob detection and converted to crane angles
`(theta1, theta2)` via `utils/pixels_to_world.m` using the bundled PnP
calibration in `data/camera_coordinates.mat`.  Boom angle `theta4` and
velocity `dtheta4` are recovered from encoder readings and smoothed with the
linear Kalman filter in `utils/linear_kalman.m`.

---

## Note on `q = 5`

The stacked trajectory sample is

```
w(i) = [ddot(theta_4), theta_1, theta_2, theta_4, dot(theta_4)]^T
```

so `q = 5`.  The first channel (`ddtheta4`) is the control input; it is
included because the behavioral framework treats inputs and outputs
symmetrically.  `theta4` and `dtheta4` are both retained so that boundary
conditions on boom position and velocity can be imposed directly in the
trajectory-generation problem.

---

## Note on `model_based/`

The `model_based/` folder contains the model-based benchmark used for
comparison in the paper.  `Verification.m` implements the 4-state load-sway
ODE driven by prescribed boom velocity and acceleration; it is called from
`MAIN.m` via `fmincon`.  This is independent of the data-driven pipeline:
the data-driven scripts use `utils/simulate_crane_one_step.m`, which embeds
the full 6-state crane dynamics with hardcoded physical parameters.

---

## Physical Parameters

These constants are hardcoded in `utils/simulate_crane_one_step.m`.

| Symbol | Value | Meaning |
|--------|-------|---------|
| `g` | 9.81 m/s² | gravity |
| `th3` | 0.8108 rad | fixed vertical boom elevation |
| `L_boom` | 2.2511 m | effective boom arm length |
| `l` | 1.0 m | rope length |
