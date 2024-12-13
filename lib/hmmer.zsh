#!/usr/bin/env zsh
# Last updated: 2024-12-03
# Tools: HMMER3
# Function libs with HMMER3
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


function search_sequences_against_Pfam_database() {
    function usage() {
        cat << EOS
Usage: search_sequences_against_Pfam_database <arg1> <arg2> <arg3>

    arg1: sp_num    <- Number of species.   (e.g., 2)
    arg2: sp_name   <- Name of species.     (e.g., "Cuscuta australis" "Cuscuta campestris")
    arg3: cpu       <- Number of CPU.       (e.g., 2)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 3 ]]; then
            usage
        fi
        sp_num=$1
        sp_names=()
        for ((i=1; i<=sp_num; i++)); do
            sp_names[i]="${2// /_}"; shift
            echo "${i}:${sp_names[i]}"
        done
        cpu=$3
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
            python3 -m fasp prefix_to_sequence_ids \
                "${taskdir}/input/${sp_name}.fasta" \
                "${taskdir}/input/${sp_name}.fasta" \
                "${sp_name}"
            cat "${taskdir}/input/${sp_name}.fasta" >> "${taskdir}/input.fasta"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        hmmscan \
            -o "${taskdir}/out.txt" \
            --domtblout "${taskdir}/domtblout.txt" \
            --cpu $cpu \
            "${DATA}/Pfam/Pfam-A.hmm" \
            "${taskdir}/input.fasta"
    }

    main "$@"
}