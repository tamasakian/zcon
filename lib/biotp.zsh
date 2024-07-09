#!/usr/bin/env zsh

#### Function libs with biotp

function construct_all_introns() {
    function usage() {
        cat << EOS
Usage:  construct_all_introns <arg1>

    arg1: genus
    arg2: org

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 2 ]]; then
            usage
        fi
        genus="$1"
        org="${2/ /_}"
    }

    function make_dir() {
        taskdir=$(make_dir_by_date $TASKFILE)
    }

    function generate_introns_gff() {
        python3 -m biotp generate_coordinate_all_introns \
            "${DATA}/${genus}/${org}.genome.gff" \
            "${taskdir}/${org}.genome.intron.gff"
    }

    function generate_introns_fasta() {
        python3 -m biotp generate_all_introns \
            "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
            "${taskdir}/${org}.genome.intron.gff" \
            "${taskdir}/${org}.intron.all.fasta"
    }

    function main() {
        parse_args "$@"
        make_dir
        generate_introns_gff
        generate_introns_fasta
    }

    main "$@"
}