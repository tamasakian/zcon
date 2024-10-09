#!/usr/bin/env zsh

#### Function libs with GS2 (Matshui et al., 2019)

function construct_tree() {
    function usage() {
        cat <<EOS
Usage: construct_tree <arg1> <arg2> <arg3> <arg4> ...

    arg1: num_orgs
    arg2: org
    arg3: num_peps
    arg4: pep
    ...

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 8 ]]; then
            usage
        fi
        
        num_orgs=$1
        shift
        orgs=()
        for ((i=1; i<=num_orgs; i++)); do
            orgs[i]="${1// /_}"
            echo "${i}:${orgs[i]}"
            shift
        done
        num_peps=$1
        shift
        typeset -g -A peps
        for ((i=1; i<=num_peps; i++)); do
            key="$1"
            shift
            value="$1"
            shift
            peps[$key]="$value"
            echo "${i}:${peps[i]}"
        done
    }

    function merge_fasta() {
        touch "${taskdir}/database.fasta"
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            cat "${DATA}/${genus}/${org}.pep.all.fasta" >> "${taskdir}/database.fasta"
        done
    }

    function generate_fasta() {
        touch "${taskdir}/pep.fasta"
        makeblastdb \
            -in "${taskdir}/database.fasta" \
            -dbtype prot \
            -hash_index \
            -parse_seqids

        for pep_key in ${(@k)peps}; do
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry "$pep_key" \
                -db "${taskdir}/database.fasta" \
                -out "$tmpfile"
            python3 -m fasp rename_header \
                "$tmpfile" \
                "$tmpfile" \
                "${peps[$pep_key]}" \
                ""
            cat "$tmpfile" >> "${taskdir}/pep.fasta"
            rm "$tmpfile"
        done
    }

    function with_gs2() {
        gs2 -e 100 -t 4 -l "${taskdir}/pep.fasta" > "${taskdir}/pep.nwk"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_fasta
        generate_fasta
        with_gs2
    }

    main "$@"
}