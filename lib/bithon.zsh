#!/usr/bin/env zsh
# Last updated: 2024-12-06
# Tools: bithon
# Function libs with bithon
: << 'FUNCTIONS'
build_longest_isoform: Build the longest isoform from NCBI/ENSEMBL genome data using bithon.
FUNCTIONS

function build_longest_isoform() {
    function usage() {
        cat << EOS
Usage: build_longest_isoform <arg1> <arg2>

    arg1: sp_num    <- Number of species.   (e.g., 2)
    arg2: sp_name   <- Name of species.     (e.g., "Cuscuta australis" "Cuscuta campestris")

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 2 ]]; then
            usage
        fi
        sp_num=$1
        sp_names=()
        echo "Species:${sp_num}"
        for ((i=1; i<=sp_num; i++)); do
            sp_names[i]="${2// /_}"; shift
            echo "${i}:${sp_names[i]}"
        done
    }

    function build_fasta() {
        mkdir "${taskdir}/input"
        for sp_name in "${sp_names[@]}"; do
            ## ENSEMBL
            if [[ -e "${DATA}/Ensembl/${sp_name}.pep.all.fasta" ]]; then
                bithon ensgls -i "${DATA}/ENSEMBL/${sp_name}.pep.all.fasta" -o "${taskdir}/input/${sp_name}.fasta" --header transcript
                continue
            fi

            ## NCBI
            mkdir "${taskdir}/input/${sp_name}"
            gn_name=${sp_name%%_*}
            cp "${DATA}/${gn_name}/${sp_name}.pep.all.fasta" "${taskdir}/input/${sp_name}/pep.fasta"
            cp "${DATA}/${gn_name}/${sp_name}.cds.all.fasta" "${taskdir}/input/${sp_name}/cds.fasta"
            bithon gls -i "${taskdir}/input/${sp_name}" -o "${taskdir}/input/${sp_name}"
            cp "${taskdir}/input/${sp_name}/longest.pep.fasta" "${taskdir}/input/${sp_name}.fasta"
            rm -r "${taskdir}/input/${sp_name}"
        done
    }

    function setup_fasta() {
        touch "${taskdir}/input.fasta"
        for sp_name in "${sp_names[@]}"; do
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                "${taskdir}/input/${sp_name}.fasta" \
                "$tmpfile" \
                "$sp_name"
            cat "$tmpfile" >> "${taskdir}/input.fasta"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        build_fasta
        setup_fasta
    }

    main "$@"
}