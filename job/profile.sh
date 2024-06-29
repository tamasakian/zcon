#!/bin/zsh

#### ZCON Profile

function declare_dir() {
    readonly ROOT=${HOME}/zcon
    readonly DATA=${ROOT}/data
    readonly DOC=${ROOT}/doc
    readonly DOWNLOAD=${ROOT}/download
    readonly JOB=${ROOT}/job
    readonly LIB=${ROOT}/lib
    readonly SCRIPT=${ROOT}/script
    readonly TASKFILE=${ROOT}/taskfile
}

function make_dir() {
    local _dirs; _dirs=($DATA $DOC $DOWNLOAD $JOB $LIB $SCRIPT $TASKFILE)
    for _dir in ${_dirs[*]}; do
        mkdir -p "$_dir"
    done
}

function declare_interpreter() {
    readonly PYTHON3=${ROOT}/venv/bin/python3.12
}

function main() {
    declare_dir
    make_dir
    declare_interpreter
}

main