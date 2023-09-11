# Two simultaneous breaks search
This folder contains the pipeline scripts required to search for two simultaneous breaks in transcripts that didn't have evidence of a single significant break.

The order to run the scripts is:
1. `prepare_transcripts.py`
2. `run_batches.py` or `run_batches_dataproc.py`
4. `merge_hts.py`
