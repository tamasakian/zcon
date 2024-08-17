#!/usr/bin/env zsh

#### Function library for SGO method

function detect_crossover_hits_by_diamond() {
    function usage() {
        cat <<EOS
Usage: detect_crossover_hits_by_diamond <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> 

    arg1: num_sgp
    arg2: org_sgp
    arg3: num_grp
    arg4: org_grp
    arg5: num_ogp
    arg6: org_ogp

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi

        num_sgp=$1
        shift
        org_sgp=()
        for ((i=1; i<=num_sgp; i++)); do
            org_sgp[i]="${1// /_}"
            echo "${i}:${org_sgp[i]}"
            shift
        done

        num_grp=$1
        shift
        org_grp=()
        for ((j=1; j<=num_grp; j++)); do
            org_grp[j]="${1// /_}"
            echo "${j}:${org_grp[j]}"
            shift
        done

        num_ogp=$1
        shift
        org_ogp=()
        for ((k=1; k<=num_ogp; k++)); do
            org_ogp[k]="${1// /_}"
            echo "${k}:${org_ogp[k]}"
            shift
        done

        total=$((num_sgp + num_grp + num_ogp))
        max_target_seqs=$((total * 10))
        echo "max_target_seqs:${max_target_seqs}"
    }

    function merge_sgp_fasta() {
        touch "${taskdir}/sgp.fasta"
        for org in "${org_sgp[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                "${DATA}/${genus}/${org}.pep.all.fasta" \
                "$tmpfile" \
                "sgp_${org}"
            cat "$tmpfile" >> "${taskdir}/sgp.fasta"
            rm "$tmpfile"
        done
    }

    function merge_grp_fasta() {
        touch "${taskdir}/grp.fasta"
        for org in "${org_grp[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                "${DATA}/${genus}/${org}.pep.all.fasta" \
                "$tmpfile" \
                "grp_${org}"
            cat "$tmpfile" >> "${taskdir}/grp.fasta"
            rm "$tmpfile"
        done
    }

    function merge_ogp_fasta() {
        touch "${taskdir}/ogp.fasta"
        for org in "${org_ogp[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                "${DATA}/${genus}/${org}.pep.all.fasta" \
                "$tmpfile" \
                "ogp"
            cat "$tmpfile" >> "${taskdir}/ogp.fasta"
            rm "$tmpfile"
        done
    }

    function makedb_reference() {
        touch "${taskdir}/reference.fasta"
        for ref in sgp grp ogp; do
            cat "${taskdir}/${ref}.fasta" >> "${taskdir}/reference.fasta"
        done
        diamond makedb \
            --in "${taskdir}/reference.fasta" \
            --db "${taskdir}/database"
    }

    function blastp_with_diamond() {
        diamond blastp \
            --db "${taskdir}/database" \
            --query "${taskdir}/rec.fasta" \
            --out "${taskdir}/hits.tsv" \
            --max-target-seqs "$max_target_seqs"
    }

    function slice_hits() {
        python3 -m biotp slice_hits_by_crossover_group \
            "${taskdir}/hits.tsv" \
            "${taskdir}/hits_slice.tsv"
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_sgp_fasta
        merge_grp_fasta
        merge_ogp_fasta
        makedb_reference
        blastp_with_diamond
        slice_hits
    }

    main "$@"
}
