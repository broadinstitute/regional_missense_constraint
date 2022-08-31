# Regional missense constraint and MPC

**WARNING: This repository is under active development!**

* Regional missense constraint (RMC): Regional variability in tolerance to missense variation within a gene
* Missense badness: Score that predicts deleteriousness of specific amino acid subsitution
* MPC (Missense badness, Polyphen-2, and Constraint): Score that predicts missense variant deleteriousness by incorporating information about the specific amino acid substitution and the region around the variant


This repository contains code to determine RMC regions in gnomAD (currently v2.1.1) and use these regions to calculate both missense badness and MPC. All methods in this repository are adapted from [Samocha et al (2017)](https://www.biorxiv.org/content/10.1101/148353v1.full.pdf).

## Repository structure:
* regional_missense_constraint/rmc/pipeline: Contains pipeline scripts that call utility functions to generate results.
* regional_missense_constraint/rmc/pipeline/two_breaks: Folder that contains scripts to search for two simultaneous breaks; see `regional_missense_constraint/rmc/pipeline/two_breaks/README.md`.
* regional_missense_constraint/rmc/utils: Contains utility functions that prepare input data or calculate results (RMC, missense badness, MPC).

## Regional missense constraint
Flowchart:
![https://lucid.app/lucidchart/f359ebaa-5c40-4fa2-9641-01b7804751a2/edit?invitationId=inv_3ae908f2-7fba-43fa-b1a2-9369a1bd8a89#]

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

## Missense badness
### Scripts
* pipeline/calculate_missense_badness.py: Pipeline code that prepares input Table and calculates missense badness score.
* utils/missense_badness.py: Utility functions that prepare data for calculations and calculates missense badness.

## MPC
### Scripts
* pipeline/calculate_mpc.py: Pipeline code that prepares input Table, calculates MPC, and annotates given input Table with MPC score.
* utils/mpc.py: Utility functions that prepare data for calculations, creates MPC model regressions, and calculates MPC.
