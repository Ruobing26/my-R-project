# Reference pipeline for the thesis

This repository contains the electronic appendix material for the bachelor thesis, including the reference pipeline, validation scripts, reference tables, and the full prompt-and-output records for all 50 model runs.

## Purpose

The repository documents how the pilot data were prepared, how the tabular prompt input was derived, and how the reference sample size tables used to benchmark the LLM outputs were computed.

The code is organised to separate the main reference pipeline, a scenario-specific validation step, and the plotting script used for the exploratory figure in Chapter 3.

## Repository contents

### Data files

- `public_D1.csv`  
  Anonymised pilot dataset D1

- `public_D2.csv`  
  Anonymised pilot dataset D2

- `public_D3.csv`  
  Anonymised pilot dataset D3

- `public_merged_sorted.csv`  
  Merged D2+D3 dataset used for model-side pooling checks

### R scripts

- `01_reference_pipeline.R`  
  Main reference pipeline. Imports the anonymised pilot datasets, reshapes them into long format, applies the complete-day restriction, computes day-specific group summaries, and calculates the reference sample sizes used to benchmark the LLM outputs.

- `02_scenario_validation_E1.R`  
  Scenario-specific validation for E1 under the exact assumptions imposed in the E1 prompt.
  
- `03_plot_prompt_input_figure.R`  
  Plotting script for the exploratory figure shown in Chapter 3.

- `04_manual_formula_checks.R`  
  Manual validation script based on the formulas presented in Chapter 2. This script reproduces selected reference calculations analytically rather than via package-based functions and compares them with the exported reference values.
  

### Output file

- `Referenztabelle_table.csv`  
  Exported reference planning table generated from the main reference pipeline

- `model_outputs/`  
  contains the full input prompt and model output for all 50 runs, organised by scenario.

- `Referenztabelle_comparison.csv`
Comparison table used for reproducibility and validation checks.

- `Referenztabelle_manual.csv`
Manual validation table for selected formula-based checks.

## Main workflow

The main reference pipeline follows these steps:

1. Import the anonymised pilot datasets
2. Reshape each dataset from wide to long format
3. Retain only complete dataset-day combinations for the prompt input and the main reference calculations
4. Compute day-specific group means, variances, and the pilot-based probabilistic effect for the Wilcoxon-Mann-Whitney approximation
5. Calculate reference sample sizes for pooled t, Welch t, pooled z, Welch z, and Wilcoxon-Mann-Whitney under one-sided and two-sided settings
6. Construct auxiliary day-to-day contrasts for broader reproducibility checks
7. Compute additional reference values for the merged D2+D3 dataset
8. Export the final reference planning table

## Notes

- For the prompt input and the main reference calculations, only complete dataset-day combinations are retained.
- The merged D2+D3 dataset is not part of the standard prompt input. It is included only to evaluate model-side aggregation decisions against a transparent reference framework.
- The plotting script assumes that the object `pilot_plus` has already been created in the main reference pipeline.
- The E1 validation script assumes that the function `get_comprehensive_n()` has already been loaded from the main reference pipeline.

## Thesis context

The repository supports the empirical analysis reported in the thesis section on data preparation, reference calculations, prompt scenarios, and output coding. It is intended to make the computational steps transparent and reproducible.
