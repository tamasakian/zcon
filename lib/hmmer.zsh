#!/usr/bin/env zsh

#### Function libs with HMMER3.4

function search_domain() {
    function usage() {
        cat << EOS
Usage: search_domain <arg1> <arg2>

    arg1: org   <- Organism name [e.g., "Arabidopsis thaliana"]
    arg2: cpu   <- Number of CPU to use [e.g., 2]

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


function search_domain_per-sequence_hits() {
    function usage() {
        cat << EOS
Usage: search_domain_per-sequence_hits <arg1>

    arg1: org   <- Organism name [e.g., "Arabidopsis thaliana"]
    arg2: cpu   <- Number of CPU to use [e.g., 2]

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


function search_domain_per-domain_hits() {
    function usage() {
        cat << EOS
Usage: search_domain_per-domain_hits <arg1>

    arg1: org   <- Organism name [e.g., "Arabidopsis thaliana"]
    arg2: cpu   <- Number of CPU to use [e.g., 2]

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
            --domtblout ${taskdir}/domtblout.txt \
            --cpu $cpu \
            ${DATA}/Pfam/Pfam-A.hmm \
            ${DATA}/${genus}/${org}.pep.all.fasta
    }

    main "$@"
}