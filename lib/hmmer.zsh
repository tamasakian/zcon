#!/usr/bin/env zsh

#### Function libs with HMMER3.4

function search_domain() {
    function usage() {
        cat << EOS
Usage: search_domain <arg1>

    arg1: org

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 1 ]]; then
            usage
        fi
        org="${1// /_}"
        genus=${org%%_*}
    }

    function search_with_hmmscan() {
        hmmscan \
            -o "${taskdir}/out.txt" \
            --tblout "${taskdir}/tblout.txt" \
            --cpu 4 \
            -E 1e-2 \
            "${DATA}/Pfam/Pfam-A.hmm" \
            "${DATA}/${genus}/${org}.pep.all.fasta"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        search_with_hmmscan
    }

    main "$@"
}