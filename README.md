# Regional missense constraint and MPC

**WARNING: This repository is under active development!**

* Regional missense constraint (RMC): Regional variability in tolerance to missense variation within a gene
* MPC (Missense deleteriousness Prediction by Constraint): RMC-based metric that predicts deleteriousness of an individual missense variant by incorporating information about the specific amino acid substitution and the region around the variant


This repository contains code to determine RMC regions in gnomAD (currently v4.1.1) and use these regions to calculate both missense badness and MPC. All methods in this repository are described in [Wang et al (2026)](https://www.biorxiv.org/content/10.1101/2024.04.11.588920v3.full.pdf).

## Repository structure:
* regional_missense_constraint/rmc/pipeline: Contains pipeline scripts that call utility functions to generate results.
* regional_missense_constraint/rmc/pipeline/two_breaks: Folder that contains scripts to search for two simultaneous breaks; see `regional_missense_constraint/rmc/pipeline/two_breaks/README.md`.
* regional_missense_constraint/rmc/utils: Contains utility functions that prepare input data or calculate results (RMC, missense badness, MPC).

## Regional missense constraint
[Flowchart](https://lucid.app/lucidchart/39bade4b-34f0-4a94-aef8-0779e36e6263/edit?viewport_loc=-506%2C-91%2C2983%2C1721%2CghqP1gcpaGqb&invitationId=inv_69b9d328-0bb2-452a-8ae6-e0ee0e5c8007)

Inputs:
* Hail Table with all possible missense variants in the human exome (called VEP context)
* gnomAD v2.1.1 public sites Hail Table

Output:
* Hail Table annotated with RMC regions and filtered to curated transcripts; removed any transcripts following this [criteria](https://gnomad.broadinstitute.org/help/why-are-constraint-metrics-missing-for-this-gene-or-annotated-with-a-note)

### Scripts:
* pipeline/regional_constraint.py: Pipeline code that prepares input Table, searches for first single/additional single breaks, and merges finalized results.
* pipeline/two_breaks: See `regional_missense_constraint/rmc/pipeline/two_breaks/README.md` for explanation of scripts in this folder.
* utils/generic.py: Utility functions that help prepare data for calculations.
* utils/constraint.py: Utility functions that calculate RMC results.

## Amino acid-level missense O/E metrics
### Notebooks
* mpc_notebooks/build_missense_badness.py: Jupyter notebook that builds missense badness, a measure of deleteriousness for an amino acid substitution class using missense O/E, normalized with nonsense and synonymous substitutions.
* mpc_notebooks/build_aa_oe_metrics: Jupyter notebook that builds the other amino acid-level missense O/E metrics described in [Wang et al (2026)](https://www.biorxiv.org/content/10.1101/2024.04.11.588920v3.full.pdf).

## MPC
### Notebooks
* mpc_notebooks/build_mpc.ipynb: Jupyter notebook that prepares the input Table to the MPC regression model, fits the MPC regression model, calculates MPC scores based on regression model fitted values, and annotates all missense variants with an MPC score.
