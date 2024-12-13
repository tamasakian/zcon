#!/usr/bin/env zsh
# Last updated: 2024-12-06
# Tools: SonicParanoid2 bithon
# Function libs with SonicParanoid2
: << 'FUNCTIONS'
identify_ortholog_groups: Identify ortholog groups.
FUNCTIONS

function identify_ortholog_groups() {
    function usage() {
        cat <<EOS
Usage: identify_ortholog_groups <arg1> <arg2> <arg3>

    arg1: sp_num    <- Number of species.        (e.g., 3)
    arg2: sp_name   <- Names  of species.        (e.g., "Cuscuta australis" "Cuscuta campestris" "Cuscuta europaea")
    arg3: threads   <- Number of threads to use. (e.g., 2)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 3 ]]; then
            usage
        fi
        sp_num=$1
        sp_names=()
        echo "Species:${sp_num}"
        for ((i=1; i<=sp_num; i++)); do
            sp_names[i]="${2// /_}"; shift
            echo "${i}:${sp_names[i]}"
        done
        threads=$2
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
        for sp_name in "${sp_names[@]}"; do
            python3 -m fasp prefix_to_sequence_ids \
                "${taskdir}/input/${sp_name}.fasta" \
                "${taskdir}/input/${sp_name}.fasta" \
                "${sp_name}"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        move_fasta_to_taskdir
        build_fasta
        setup_fasta
        mkdir "${taskdir}/output"
        sonicparanoid -i "${taskdir}/input" -o "${taskdir}/output" -t $threads
    }

    main "$@"
}