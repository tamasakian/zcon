#!/usr/bin/env zsh

#### Function libs with OrthoFinder

function search_orthogroups() {
    function usage() {
        cat <<EOS
Usage search_orthogroups <arg1> <arg2> ...

    arg1: num (more than 3)
    arg2: org
    ...

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 4 ]]; then
            usage
        fi
        
        num=$1
        shift
        orgs=()
        for ((i=1; i<=num; i++)); do
            orgs[i]="${1// /_}"
            echo "${i}:${orgs[i]}"
            shift
        done
    }

    function move_fasta_to_taskdir() {
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cp "${DATA}/${genus}/${org}.pep.all.fasta" "${taskdir}/${org}.fasta"
        done
    }

    function search_orthogroups_with_orgs() {
        oethofinder \
            -t 128 \
            -a 16 \
            -f "$taskdir" \
            -T "raxml-ng"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        move_fasta_to_taskdir
        search_orthogroups_with_orgs
    }

    main "$@"
}