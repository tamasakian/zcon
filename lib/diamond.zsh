#!/usr/bin/env zsh

#### Function libs with DIAMOND

function output_blastp_by_rgo() {
    ## rgo is;
    ## r: recipient, rec;
    ## g: group, grp;
    ## o: outgroup, ogp
    function usage() {
        cat <<EOS
Usage: output_blastp_by_3_genus <arg1> (<arg2> ... <arg(1+x)>) <arg(2+x)> (<arg(3+x)> ... <arg(2+x+y)>) <arg(3+x+y)> (<arg(4+x+y)> ... <arg(3+x+y+z)>) <arg(4+x+y+z)> 

    arg1: num_rec [x]
    arg2: org_rec
    ...

    arg(2+x): num_grp [y]
    arg(3+x): org_grp
    ...

    arg(3+x+y): num_ogp [z]
    arg(4+x+y): org_ogp
    ...

    arg(4+x+y+z): evalue

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi

        num_rec=$1
        shift
        org_recs=()
        for ((i=1; i<=num_rec; i++)); do
            org_recs[i]="${1// /_}"
            echo "${i}:${org_recs[i]}"
            shift
        done

        num_grp=$1
        shift
        org_grps=()
        for ((j=1; j<=num_grp; j++)); do
            org_grps[j]="${1// /_}"
            echo "${j}:${org_grps[j]}"
            shift
        done

        num_ogp=$1
        shift
        org_ogps=()
        for ((k=1; k<=num_ogp; k++)); do
            org_ogps[k]="${1// /_}"
            echo "${k}:${org_ogps[k]}"
            shift
        done

        evalue=$1
        echo "evalue:${evalue}"
    }

    function merge_rec_fasta() {
        touch "${taskdir}/rec.fasta"
        for org in "${org_recs[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m biotp prefix_to_headers \
                "${DATA}/${genus}/${org}.pep.all.fasta" \
                "$tmpfile" \
                "rec"
            cat "$tmpfile" >> "${taskdir}/rec.fasta"
            rm "$tmpfile"
        done
    }

    function merge_grp_fasta() {
        touch "${taskdir}/grp.fasta"
        for org in "${org_grps[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m biotp prefix_to_headers \
                "${DATA}/${genus}/${org}.pep.all.fasta" \
                "$tmpfile" \
                "grp"
            cat "$tmpfile" >> "${taskdir}/grp.fasta"
            rm "$tmpfile"
        done
    }

    function merge_ogp_fasta() {
        touch "${taskdir}/ogp.fasta"
        for org in "${org_ogps[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m biotp prefix_to_headers \
                "${DATA}/${genus}/${org}.pep.all.fasta" \
                "$tmpfile" \
                "ogp"
            cat "$tmpfile" >> "${taskdir}/ogp.fasta"
            rm "$tmpfile"
        done
    }

    function makedb_reference() {
        touch "${taskdir}/reference.fasta"
        for ref in rec grp ogp; do
            cat "${ref}.fasta" >> "${taskdir}/reference.fasta"
        done
        diamond makedb \
            --in "${taskdir}/reference.fasta" \
            -d "${taskdir}/database"
    }

    function blastp_with_diamond() {
        touch matches.tsv
        diamond blastp \
            -d "${taskdir}/database" \
            -q "${taskdir}/rec.fasta" \
            -o "${taskdir}/matches.tsv"
    }



    function main() {
        parse_args "$@"
        make_taskdir
        merge_rec_fasta
        merge_grp_fasta
        merge_ogp_fasta
        makedb_reference
        blastp_with_diamond
    }

    main "$@"
}