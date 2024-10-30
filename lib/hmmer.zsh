#!/usr/bin/env zsh
# Last updated: 2024-10-30

#### Function libs with HMMER3
: << 'FUNCTIONS'
search_domain:              Search for protein domains using hmmscan.
search_domain_per_sequence: Search for protein domains using hmmscan and Output domain hits per-sequence.
search_domain_per_domain:   Search for protein domains using hmmscan and Output domain hits per-domain.
FUNCTIONS


function search_domain() {
    function usage() {
        cat << EOS
Usage: search_domain <arg1> <arg2>

    arg1: org   <- Organism name.        (e.g., "Arabidopsis thaliana")
    arg2: cpu   <- Number of CPU to use. (e.g., 2)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 2 ]]; then
            usage
        fi
        org="${1// /_}"
        cpu=$2
        genus=${org%%_*}
    }

    function main() {
        parse_args "$@"
        make_taskdir
        hmmscan \
            -o ${taskdir}/out.txt \
            --cpu $cpu \
            ${DATA}/Pfam/Pfam-A.hmm \
            ${DATA}/${genus}/${org}.pep.all.fasta
    }

    main "$@"
}


function search_domain_per_sequence() {
    function usage() {
        cat << EOS
Usage: search_domain_per_sequence <arg1> <arg2>

    arg1: org   <- Organism name.        (e.g., "Arabidopsis thaliana")
    arg2: cpu   <- Number of CPU to use. (e.g., 2)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 2 ]]; then
            usage
        fi
        org="${1// /_}"
        cpu=$2
        genus=${org%%_*}
    }

    function main() {
        parse_args "$@"
        make_taskdir
        hmmscan \
            -o ${taskdir}/out.txt \
            --tblout ${taskdir}/tblout.txt \
            --cpu $cpu \
            ${DATA}/Pfam/Pfam-A.hmm \
            ${DATA}/${genus}/${org}.pep.all.fasta
    }

    main "$@"
}


function search_domain_per_domain() {
    function usage() {
        cat << EOS
Usage: search_domain_per_domain <arg1> <arg2>

    arg1: org   <- Organism name.        (e.g., "Arabidopsis thaliana")
    arg2: cpu   <- Number of CPU to use. (e.g., 2)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 2 ]]; then
            usage
        fi
        org="${1// /_}"
        cpu=$2
        genus=${org%%_*}
    }

    function main() {
        parse_args "$@"
        make_taskdir
        hmmscan \
            -o ${taskdir}/out.txt \
            --domtblout ${taskdir}/domtblout.txt \
            --cpu $cpu \
            ${DATA}/Pfam/Pfam-A.hmm \
            ${DATA}/${genus}/${org}.pep.all.fasta
    }

    main "$@"
}