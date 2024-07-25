#!/usr/bin/env zsh

#### Function libs with DIAMOND

function output_blastp_by_rgo() {
    function usage() {
        cat <<EOS
Usage: output_blastp_by_rgo <arg1> (<arg2> ... <arg(1+x)>) <arg(2+x)> (<arg(3+x)> ... <arg(2+x+y)>) <arg(3+x+y)> (<arg(4+x+y)> ... <arg(3+x+y+z)>) <arg(4+x+y+z)> 

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
            --out "${taskdir}/matches.tsv" \
            --max-target-seqs 25
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

function search_hgt_by_rgo() {
    ## This function is based on output_blastp_by_rgo
    function usage() {
        cat <<EOS
Usage: search_hgt_by_rgo <arg1> (<arg2> ... <arg(1+x)>) <arg(2+x)> (<arg(3+x)> ... <arg(2+x+y)>) <arg(3+x+y)> (<arg(4+x+y)> ... <arg(3+x+y+z)>) <arg(4+x+y+z)> 

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

    function slice_hgt() {
        python3 -m biotp slice_rgo_by_hgt \
            "${taskdir}/matches.tsv" \
            "${taskdir}/hgt_matches.tsv" \
            60
    }

    function main() {
        output_blastp_by_rgo "$@"
        slice_hgt
    }

    main "$@"
}

function detect_hgt_by_rgo() {
    ## This function is based on search_hgt_by_rgo
    function usage() {
        cat <<EOS
Usage: detect_hgt_by_rgo <arg1> (<arg2> ... <arg(1+x)>) <arg(2+x)> (<arg(3+x)> ... <arg(2+x+y)>) <arg(3+x+y)> (<arg(4+x+y)> ... <arg(3+x+y+z)>) <arg(4+x+y+z)> 

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

    function merge_rec_cds_fasta() {
        touch "${taskdir}/rec.cds.fasta"
        for org in "${org_recs[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m biotp rename_headers_feature \
                "${DATA}/${genus}/${org}.cds.all.fasta" \
                "$tmpfile1" \
                "protein_id"
            python3 -m biotp prefix_to_headers \
                "$tmpfile1" \
                "$tmpfile2" \
                "rec"
            cat "$tmpfile2" >> "${taskdir}/rec.cds.fasta"
            rm "$tmpfile1" "$tmpfile2"
        done
    }

    function merge_grp_cds_fasta() {
        touch "${taskdir}/grp.cds.fasta"
        for org in "${org_grps[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m biotp rename_headers_feature \
                "${DATA}/${genus}/${org}.cds.all.fasta" \
                "$tmpfile1" \
                "protein_id"
            python3 -m biotp prefix_to_headers \
                "$tmpfile1" \
                "$tmpfile2" \
                "grp"
            cat "$tmpfile2" >> "${taskdir}/grp.cds.fasta"
            rm "$tmpfile1" "$tmpfile2"
        done
    }

    function merge_ogp_cds_fasta() {
        touch "${taskdir}/ogp.cds.fasta"
        for org in "${org_ogps[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m biotp rename_headers_feature \
                "${DATA}/${genus}/${org}.cds.all.fasta" \
                "$tmpfile1" \
                "protein_id"
            python3 -m biotp prefix_to_headers \
                "$tmpfile1" \
                "$tmpfile2" \
                "ogp"
            cat "$tmpfile2" >> "${taskdir}/ogp.cds.fasta"
            rm "$tmpfile1" "$tmpfile2"
        done
    }

    function merge_reference() {
        touch "${taskdir}/reference.cds.fasta"
        for ref in rec grp ogp; do
            cat "${taskdir}/${ref}.cds.fasta" >> "${taskdir}/reference.cds.fasta"
        done
    }

    function slice_fasta() {
        _recs=($(cut -f 1 "${taskdir}/hgt_matches.tsv" | sort -u))
        python3 -m biotp slice_records_by_names \
            "${taskdir}/rec.cds.fasta" \
            "${taskdir}/qry.cds.fasta" \
            "${_recs[@]}"
        _references=($(cut -f 2 "${taskdir}/hgt_matches.tsv" | sort -u))
        python3 -m biotp slice_records_by_names \
            "${taskdir}/reference.cds.fasta" \
            "${taskdir}/ref.cds.fasta" \
            "${_references[@]}"
    }

    function makeblastdb_reference() {
        makeblastdb \
            -in "${taskdir}/ref.cds.fasta" \
            -dbtype nucl -hash_index -parse_seqids
    }

    function detect_hgt_by_blastn() {
        blastn \
            -outfmt 6 -evalue 10 \
            -db "${taskdir}/ref.cds.fasta" \
            -query "${taskdir}/qry.cds.fasta" \
            -out "${taskdir}/results.tsv"
        python3 -m biotp detect_hgt \
            "${taskdir}/results.tsv" \
            "${taskdir}/hgt.tsv"
    }

    function main() {
        search_hgt_by_rgo "$@"
        merge_rec_cds_fasta
        merge_grp_cds_fasta
        merge_ogp_cds_fasta
        merge_reference
        slice_fasta
        makeblastdb_reference
        detect_hgt_by_blastn
    }

    main "$@"
}