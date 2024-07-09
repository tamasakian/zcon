#!/usr/bin/env zsh

#### ZCON Profile

export ROOT=${HOME}/zcon
export DATA=${ROOT}/data
export DOC=${ROOT}/doc
export DOWNLOAD=${ROOT}/download
export JOB=${ROOT}/job
export LIB=${ROOT}/lib
export SCRIPT=${ROOT}/script
export TASKFILE=${ROOT}/taskfile

for _dir in $DATA $DOC $DOWNLOAD $JOB $LIB $SCRIPT $TASKFILE; do
    mkdir -p "$_dir"
done

source "${ROOT}/venv/bin/activate"