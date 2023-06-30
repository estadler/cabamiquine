# `analysis/` directory

holds scripts used to transform analyse clean data and create outputs stored in `output/`

Files:

- `01_estimates-from-data.R` estimates the frequency of resistant mutants in the data
- `02_data-estimates-comparison.R` compares the frequency estimates of resistant parasites from different data
- `03_deterministic-model.R` runs the deterministic model for the frequency of resistant mutants by generation
- `04_stochastic-model.R` runs simulations of the stochastic model for the frequency of resistant mutants
- `05_stochastic-sim-experiments.R` runs stochastic simulations of the different experiments with experiment-specific parameter values and estimates the frequency of resistant mutants as in the data
- `06_visualisations.R` produces figures and plots of the data and simulation estimates of the frequency of resistant parasites, the sensitivity analysis of the estimates from the simulated experiments, and the probability of emergence of resistant parasites only after treatment