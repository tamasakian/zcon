#!/usr/bin/env zsh

#### Function libs

function make_taskdir() {
    ## This function makes a task directory.
    function main() {
        dirname=$(date +"%Y-%m-%d-%H-%M-%S")
        taskdir=${TASKFILE}/${dirname}
        mkdir -p $taskdir
    }

    main
}