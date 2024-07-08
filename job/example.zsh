#!/usr/bin/env zsh

#### ZCON Example

source $(dirname "$0")/profile.zsh

## Library
function import_libs() {
    source ${LIB}/function.zsh
    source ${LIB}/datasets.zsh
    source ${LIB}/blast.zsh
    source ${LIB}/mcscan.zsh
    source ${LIB}/mummer.zsh
    source ${LIB}/meme.zsh
}

## Job
function job_19221111() {
    echo "Happy Birthday!"
    fetch_genome_by_genus "Homo"
}

## Run
function main() {
    import_libs
    job_19221111
}

main