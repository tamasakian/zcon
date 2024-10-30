#!/usr/bin/env zsh
# Last updated: 2024-10-30

#### Function libs with SonicParanoid2 (Cosentino et al., 2024)
: << 'FUNCTIONS'
identify_ortholog_groups: Identify ortholog groups.
FUNCTIONS

function identify_ortholog_groups() {
    function usage() {
        cat <<EOS
Usage: identify_ortholog_groups <arg1> <arg2> <arg3>

    arg1: num       <- Number of organisms.      (e.g., 4)
    arg2: org       <- Names  of organisms.      (e.g., "Cuscuta australis" "Cuscuta campestris" ... "Cuscuta europaea")
    arg3: threads   <- Number of threads to use. (e.g., 2)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 5 ]]; then
            usage
        fi
        num=$1; shift
        orgs=()
        for ((i=1; i<=num; i++)); do
            orgs[i]="${1// /_}"; shift
        done
        threads=$1; shift
    }

    function move_fasta_to_taskdir() {
        mkdir ${taskdir}/input ${taskdir}/output
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m fasp exclude_isoforms_by_length \
                ${DATA}/${genus}/${org}.pep.all.fasta \
                $tmpfile1 \
                ${DATA}/${genus}/${org}.genome.gff 
            python3 -m fasp exclude_non_nuclear_proteins \
                $tmpfile1 \
                $tmpfile2
            python3 -m fasp prefix_to_sequence_ids \
                $tmpfile2 \
                ${taskdir}/input/${org}.fasta \
                ${org}
            rm $tmpfile1 $tmpfile2
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        move_fasta_to_taskdir
        sonicparanoid -i ${taskdir}/input -o ${taskdir}/output -t $threads
    }

    main "$@"
}