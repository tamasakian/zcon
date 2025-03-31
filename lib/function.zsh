#!/usr/bin/env zsh
# Last updated: 2024-10-30

#### Function libs
: << 'FUNCTIONS'
make_taskdir: Creates a task directory with a timestamp.
FUNCTIONS


function make_taskdir() {
    dirname=$(date +"%Y-%m-%d-%H-%M-%S")
    taskdir=${TASKFILE}/${dirname}
    mkdir -p $taskdir
}