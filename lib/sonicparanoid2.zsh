#!/usr/bin/env zsh

#### Function libs with SonicParanoid2 (Salvatore et al., 2024)

function infer_orthology() {
    function usage() {
        cat <<EOS
Usage: infer_orthology <arg1> <arg2> ...

    arg1: num
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
        mkdir "${taskdir}/input" "${taskdir}/output"
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
                "${taskdir}/input/${org}.fasta" \
                "${org}_"
            rm "$tmpfile1" "$tmpfile2"
        done
    }

    function sonicparanoid_with_apptainer() {
        module load apptainer
        apptainer exec /usr/local/biotools/s/sonicparanoid:2.0.8--py312h1f1cfbb_0 \
            sonicparanoid -i "${taskdir}/input" -o "${taskdir}/output" -t 16
    }

    function main() {
        parse_args "$@"
        make_taskdir
        move_fasta_to_taskdir
        sonicparanoid_with_apptainer
    }

    main "$@"
}