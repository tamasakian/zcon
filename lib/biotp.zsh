#!/usr/bin/env zsh
# Last updated: 2024-12-04
# Tools: BIOTP FASP
# Function libs with BIOTP
: << 'FUNCTIONS'
generate_all_introns_fasta: Generate an intron FASTA file. 
calculate_all_introns_lengths: Calculate lengths from the generated intron FASTA file. 
FUNCTIONS

function generate_all_introns_fasta() {
    function usage() {
        cat << EOS
Usage: generate_all_introns_fasta <arg1>

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

    function generate_gff() {
        python3 -m biotp generate_introns \
            "${DATA}/${genus}/${org}.genome.gff" \
            "${taskdir}/${org}.intron.gff"
    }

    function generate_fasta() {
        python3 -m fasp generate_introns \
            "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
            "${taskdir}/${org}.intron.all.fasta" \
            "${taskdir}/${org}.intron.gff"
    }   

    function main() {
        parse_args "$@"
        make_taskdir
        generate_gff
        generate_fasta
    }

    main "$@"
}

function calculate_all_introns_lengths() {
    function usage_2() {
        cat << EOS
Usage: calculate_all_introns_lengths <arg1>

    arg1: org

EOS
        exit 1
    }

    function parse_args_2() {
        if [[ $# != 1 ]]; then
            usage_2
        fi
        org="${1// /_}"
        genus=${org%%_*}
    }

    function calculate_lengths() {
        python3 -m fasp measure_lengths \
            "${taskdir}/${org}.intron.all.fasta" \
            "${taskdir}/${org}.intron.length.tsv"
    }

    function main() {
        generate_all_introns_fasta "$@"
        calculate_lengths
    }

    main "$@"
}

function make_dict_pepid_by_gff() {
    function usage() {
        cat << EOS
Usage: make_dict_pepid_by_gff <arg1> <arg2> <arg3> <arg4> <arg5>
    
    arg1: genus
    arg2: org
    arg3: kind
    arg4: key
    arg5: pattern

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 5 ]]; then
            usage
        fi
        genus="$1"
        org="${2/ /_}"
        kind="$3"
        key="$4"
        pattern="$5"
    }

    function make_dict() {
        python3 -m biotp make_dict_pepid \
            "${DATA}/${genus}/${org}.genome.gff" \
            "${taskdir}/${org}.dict_pepid.csv" \
            "$kind" \
            "$key" \
            "$pattern"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        make_dict
    }

    main "$@"

}

function generate_all_upstream_reagions_fasta() {
    function usage() {
        cat << EOS
Usage: generate_all_upstream_reagions_fasta <arg1> <arg2>

    arg1: org
    arg2: bp

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 2 ]]; then
            usage
        fi
        org="${1// /_}"
        bp=$2
        genus=${org%%_*}
    }

    function generate_gff() {
        python3 -m biotp generate_upstream_regions \
            "${DATA}/${genus}/${org}.genome.gff" \
            "${taskdir}/${org}.upstream_${bp}.gff" \
            "$bp"
    }

    function generate_fasta() {
        python3 -m fasp generate_upstream_regions \
            "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
            "${taskdir}/${org}.upstream_${bp}.all.fasta" \
            "${taskdir}/${org}.upstream_${bp}.gff"
    }   

    function main() {
        parse_args "$@"
        make_taskdir
        generate_gff
        generate_fasta
    }

    main "$@"
}


function generate_all_downstream_reagions_fasta() {
    function usage() {
        cat << EOS
Usage: generate_all_downstream_reagions_fasta <arg1> <arg2>

    arg1: org
    arg2: bp

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 2 ]]; then
            usage
        fi
        org="${1// /_}"
        bp=$2
        genus=${org%%_*}
    }

    function generate_gff() {
        python3 -m biotp generate_downstream_regions \
            "${DATA}/${genus}/${org}.genome.gff" \
            "${taskdir}/${org}.downstream_${bp}.gff" \
            "$bp"
    }

    function generate_fasta() {
        python3 -m fasp generate_downstream_regions \
            "${DATA}/${genus}/${org}.dna.toplevel.fasta" \
            "${taskdir}/${org}.downstream_${bp}.all.fasta" \
            "${taskdir}/${org}.downstream_${bp}.gff"
    }   

    function main() {
        parse_args "$@"
        make_taskdir
        generate_gff
        generate_fasta
    }

    main "$@"
}