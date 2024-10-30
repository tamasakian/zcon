#!/usr/bin/env zsh
# Last updated: 2024-10-30

#### Function libs
# make_taskir   <- Creates a task directory with a timestamp.

function make_taskdir() {
    dirname=$(date +"%Y-%m-%d-%H-%M-%S")
    taskdir=${TASKFILE}/${dirname}
    mkdir -p $taskdir
}