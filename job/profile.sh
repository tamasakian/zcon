#!/usr/bin/env zsh

setopt errexit
setopt nounset
setopt pipefail

#### ZCON Profile

function declare_dirs() {
    export ROOT=${HOME}/zcon
    export DATA=${ROOT}/data
    export DOC=${ROOT}/doc
    export DOWNLOAD=${ROOT}/download
    export JOB=${ROOT}/job
    export LIB=${ROOT}/lib
    export SCRIPT=${ROOT}/script
    export TASKFILE=${ROOT}/taskfile
}

function make_dirs() {
    local _dirs; _dirs=($DATA $DOC $DOWNLOAD $JOB $LIB $SCRIPT $TASKFILE)
    for _dir in ${_dirs[*]}; do
        mkdir -p "$_dir"
    done
}

function main() {
    declare_dirs
    make_dirs
}

main

alias pyvenv='source ${ROOT}/venv/bin/activate && ${ROOT}/venv/bin/python3'
