#!/usr/bin/env zsh

#### ZCON Example

source $(dirname "$0")/profile.sh

## Library
function import_lib() {
    source ${LIB}/function.sh
    source ${LIB}/datasets.sh
    source ${LIB}/blast.sh
    source ${LIB}/mcscan.sh
    source ${LIB}/mummer.sh
    source ${LIB}/meme.sh
}

## Job
function job_19221111() {
    echo "Happy Birthday!"
    fetch_genome_by_genus "Homo"
}

## Run
function main() {
    import_lib
    job_19221111
}

main