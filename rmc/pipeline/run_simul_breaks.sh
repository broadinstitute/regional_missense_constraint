#!/bin/bash
usage() {
cat << EOF
    This script runs `simultaneous_breaks.py` on groups of transcripts serially.

    Inputs:
        -b    Google cloud bucket containing TSVs with groups of transcripts.
        -p    Prefix for transcript TSV file name.
        -c    Google cloud cluster name.
        -t    Chi square threshold.
        -w    Minimum window size to search for simultaneous breaks.
        -s    Path to Python script to run.

    Example command:
    bash pipeline/run_simul_breaks.sh -b gs://regional_missense_constraint/temp/transcripts -p xa -c kc -t 9.2 -w 100 -s /Users/kchao/code/regional_missense_constraint/rmc/utils/simultaneous_breaks.py

EOF
}

# Check num of command line args
# NOTE: This number is 12 here because it counts both the flag and the value
# e.g., counts both the flag "-w" and the value "100"
if [[ $# -lt 12 ]]; then
    usage
    exit 0
fi

# Parse command line args
while getopts "b:p:c:t:w:s:h" opt; do
    case $opt in
        b)
            bucket=$OPTARG
        ;;
        p)
            prefix=$OPTARG
        ;;
        c)
            cluster=$OPTARG
        ;;
        t)
            chisq=$OPTARG
        ;;
        w)
            window_size=$OPTARG
        ;;
        s)
            python_path=$OPTARG
        ;;
        h)
            usage
            exit 0
        ;;
        \?)
            usage
            exit 0
        ;;
    esac
done


# Cycle through lists of transcript TSVs
transcript_tsv_list=$(gsutil ls ${bucket}/${prefix}*)
for tsv in ${transcript_tsv_list[@]}; do
    tsv_name=$(basename ${tsv})
    log=${tsv_name}.log
    echo "Starting job for transcripts in ${tsv_name}. Log file: ${log}"
    (
    hailctl dataproc submit ${cluster} ${python_path} --min-window-size ${window_size} --chisq-threshold ${chisq} --transcript-tsv ${tsv}
    wait
    ) &> $log
done
