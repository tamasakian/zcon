#!/usr/bin/env zsh

#### Function libs with OrthoFinder

function search_orthogroups() {
    function usage() {
        cat <<EOS
Usage: search_orthogroups <arg1> <arg2> ...

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
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            genus=${org%%_*}
            python3 -m fasp exclude_isoforms_by_length \
                "${DATA}/${genus}/${org}.pep.all.fasta" \
                "$tmpfile1" \
                "${DATA}/${genus}/${org}.genome.gff" 
            python3 -m fasp exclude_non_nuclear_proteins \
                "$tmpfile1" \
                "$tmpfile2"
            python3 -m fasp prefix_to_sequence_ids \
                "$tmpfile2" \
                "${taskdir}/${org}.fasta" \
                "${org}"
            rm "$tmpfile1" "$tmpfile2"
        done
    }

    function search_orthogroups_with_orgs() {
        orthofinder -f "$taskdir"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        move_fasta_to_taskdir
        search_orthogroups_with_orgs
    }

    main "$@"
}