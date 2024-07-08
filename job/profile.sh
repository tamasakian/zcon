#!/usr/bin/env zsh

#### ZCON Profile

setopt err_exit
setopt no_unset
setopt pipe_fail

export ROOT=${HOME}/zcon
export DATA=${ROOT}/data
export DOC=${ROOT}/doc
export DOWNLOAD=${ROOT}/download
export JOB=${ROOT}/job
export LIB=${ROOT}/lib
export SCRIPT=${ROOT}/script
export TASKFILE=${ROOT}/taskfile

local _dirs; _dirs=($DATA $DOC $DOWNLOAD $JOB $LIB $SCRIPT $TASKFILE)
for _dir in ${_dirs[*]}; do
    mkdir -p "$_dir"
done

alias pyvenv='source ${ROOT}/venv/bin/activate && ${ROOT}/venv/bin/python3'
