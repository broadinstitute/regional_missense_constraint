# Two simultaneous breaks search
This folder contains the pipeline scripts required to search for two simultaneous breaks in transcripts that didn't have evidence of a single significant break.

The order to run the scripts is:
1. `prepare_transcripts.py`
2. `run_batches_dataproc.py`
4. `merge_hts.py`

The search also ran in Hail Batch for gnomAD v2, through a `run_batches.py` that carried
its own copy of the search functions because Batch PythonJobs can't import this repo.
That script is gone; see the `rmc_gnomad_v2` branch or git history if the Batch path is
ever needed again.
