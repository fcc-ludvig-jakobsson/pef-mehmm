# Mixed-Effects Hidden Markov Model for PEF (SAEM Implementation)
This repository contains an R implementation of a mixed-effects hidden Markov model (MEHMM) for modeling home-measured peak expiratory flow (PEF), together with a full simulation and estimation pipeline based on the SAEM (Stochastic Approximation EM) algorithm.

The model is designed to extract interpretable features of lung function dynamics and quantify treatment effects from longitudinal respiratory data.

## Background
Home-measured peak expiratory flow (PEF) provides dense time series data reflecting lung function. These data often exhibit:
- High variability
- Sustained periods of worsened lung function
- Individual-specific dynamics.

Standard approaches may not fully capture this temporal structure. This project instead uses a mixed-effects hidden Markov model, where:
- Observations are generated from latent (hidden) disease states
- Each subject has individual-specific parameters (random effects)
- Treatment effects can influence all model parameters

The project resulted in an article published in CPT:PSP
- DOI: 10.1002/psp4.70281

## Hidden Markov Model
The model assumes two hidden states:
- State 0: High PEF
- State 1: Low PEF

At each time point:

- Observations follow a normal distribution conditioned on state
- Transitions between states follow a Markov process

### Parameters
Each individual has 5 parameters:

- Baseline PEF (μ₀)
- Relative drop between states (d)
- Within-state variance (σ²)
- Transition probability (0 → 1)
- Transition probability (1 → 0)

### Mixed Effects

Parameters vary across individuals via random effects

Population parameters: μ, β (treatment), ω (variance)

### Dose-Response Models
Supported treatment effects include:

- none
- const (binary effect)
- categorical
- linear
- emax

## Key Functions
### Simulation

- `simulate_two_state_model()` \
  Simulates individual time series with hidden states

### Inference
- `SAEM()` \
Main estimation routine (population parameters)

- `MCMC()` \
Individual-level sampling (within SAEM)


- `forward_algorithm()` \
Computes likelihood via HMM forward recursion


### Post-processing

- `viterbi()` → Most likely hidden state sequence
- `MAP()` → Posterior mean estimates (approximate MAP)
- `IMP_log_likelihood()` → Importance sampling likelihood
- `info_matrix_estimation()` → Fisher information → standard errors
- `calc_confusion_matrix()` → State classification accuracy

## Correspondance
Please contact ludvig.jakobsson@fcc.chalmers.se with questions.
