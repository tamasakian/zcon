#!/usr/bin/env zsh

#### Function libs

function make_dir_by_date() {
    function usage() {
        cat <<EOS
Usage:  make_dir_by_date <arg1>

    arg1: dir

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 1 ]]; then
            usage
        fi
        dir="$1"
    }

    function make_dir() {
        dirname=$(date +"%Y-%m-%d-%H-%M-%S")
        mkdir -p "${dir}/${dirname}"
    }

    function output_dirname() {
        echo "${dir}/${dirname}"
    }

    function main() {
        parse_args "$@"
        make_dir
        output_dirname
    }

    main "$@"
}

function make_taskdir() {
    ## This function is based on make_dir_by_date
    function usage() {
        cat <<EOS
Usage: make_taskdir

    This fucntion does not need any args.

EOS
        exit 1
    }

    function main() {
        taskdir=$(make_dir_by_date $TASKFILE)
    }

    main
}