#!/usr/bin/env bash

# manually add your path if you're on Windows (don't forget escape character on \, \\)
START_DIR=$(pwd)
RESOURCE_DIR=$START_DIR
FILE_OUT_DIR=$START_DIR
FILE_NAME="RenyiPhaseImpDecoy46Results";


TOTAL_JOBS=16
TOTAL_PARALLEL_JOBS=16

if [[ "$FILE_OUT_DIR" == /* ]]; then
    LOG_DIR=$FILE_OUT_DIR # abs
else
    LOG_DIR="${START_DIR}/${FILE_OUT_DIR}" # rel
fi

echo "Start Directory:    " $START_DIR
echo "Resource Directory: " $RESOURCE_DIR
echo "Save Directory:     " $FILE_OUT_DIR
echo

# note every variable statement except {JOB_IND} is replaced before xargs
# replaces {JOB_IND}.
seq 1 "$TOTAL_JOBS" | xargs -P "$TOTAL_PARALLEL_JOBS" -I {JOB_IND} bash -c "
    echo \"starting Job: {JOB_IND}\"
    
    # You may have to remove '-nojvm' if you encounter esoteric error messages.
    matlab -noFigureWindows -singleCompThread -sd \"${START_DIR}\" \
        -batch \"addpath(genpath('${RESOURCE_DIR}')); \
        RenyiPhaseImpDecoy46batch('${FILE_OUT_DIR}', '${FILE_NAME}', \
        {JOB_IND}, ${TOTAL_JOBS});\" > \"${LOG_DIR}/${FILE_NAME}_{JOB_IND}.log\" \
        2>&1
"

echo "Jobs Finished!"