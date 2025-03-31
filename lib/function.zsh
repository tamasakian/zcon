#!/usr/bin/env zsh
# Last updated: 2025-03-31

#### Function libs
: << 'FUNCTIONS'
make_taskdir: Creates a task directory with a timestamp.
FUNCTIONS


function make_taskdir() {
    dirname=$(date +"%Y-%m-%d-%H-%M-%S")
    taskdir=${TASKFILE}/${dirname}
    mkdir -p $taskdir
}