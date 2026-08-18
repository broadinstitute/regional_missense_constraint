# Two simultaneous breaks search
This folder contains the pipeline scripts required to search for two simultaneous breaks in transcripts that didn't have evidence of a single significant break.

The order to run the scripts is:
1. `prepare_transcripts.py`
2. `run_batches_dataproc.py`
4. `merge_hts.py`

`run_batches.py` runs the search in Hail Batch and was only used for the gnomAD v2
searches. It duplicates the search functions from `rmc/utils/simultaneous_breaks.py`
because Batch PythonJobs can't import this repo, and the two copies have since drifted.
