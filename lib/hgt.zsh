#!/usr/bin/env zsh

#### Function library for HGT research

function detect_crossover_besthits() {
    function usage() {
        cat <<EOS
Usage: detect_crossover_besthits <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> 

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
                "sgp_${org}_"
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
                "grp_${org}_"
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
                "ogp_${org}_"
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
            --query "${taskdir}/sgp.fasta" \
            --out "${taskdir}/hits.tsv" \
            --max-target-seqs "$max_target_seqs"
    }

    function slice_hits() {
        python3 -m biotp slice_hits_by_crossover_group \
            "${taskdir}/hits.tsv" \
            "${taskdir}/hits_slice.tsv"
        python3 -m biotp output_besthit_for_subgroup \
            "${taskdir}/hits.tsv" \
            "${taskdir}/besthits.tsv"
    }

    function merge_sgp_fasta_2() {
        touch "${taskdir}/sgp.cds.all.fasta"
        for org in "${org_sgp[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m fasp rename_headers_feature \
                "${DATA}/${genus}/${org}.cds.all.fasta" \
                "$tmpfile1" \
                "protein_id"
            python3 -m fasp prefix_to_sequence_ids \
                "$tmpfile1" \
                "$tmpfile2" \
                "sgp_${org}_"
            cat "$tmpfile2" >> "${taskdir}/sgp.cds.all.fasta"
            rm "$tmpfile1" "$tmpfile2"
        done
    }

    function merge_grp_fasta_2() {
        touch "${taskdir}/grp.cds.all.fasta"
        for org in "${org_grp[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m fasp rename_headers_feature \
                "${DATA}/${genus}/${org}.cds.all.fasta" \
                "$tmpfile1" \
                "protein_id"
            python3 -m fasp prefix_to_sequence_ids \
                "$tmpfile1" \
                "$tmpfile2" \
                "grp_${org}_"
            cat "$tmpfile2" >> "${taskdir}/grp.cds.all.fasta"
            rm "$tmpfile1" "$tmpfile2"
        done
    }

    function merge_ogp_fasta_2() {
        touch "${taskdir}/ogp.cds.all.fasta"
        for org in "${org_ogp[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m fasp rename_headers_feature \
                "${DATA}/${genus}/${org}.cds.all.fasta" \
                "$tmpfile1" \
                "protein_id"
            python3 -m fasp prefix_to_sequence_ids \
                "$tmpfile1" \
                "$tmpfile2" \
                "ogp_${org}_"
            cat "$tmpfile2" >>"${taskdir}/ogp.cds.all.fasta"
            rm "$tmpfile1" "$tmpfile2"
        done
    }

    function makeblastdb_reference() {
        touch "${taskdir}/reference.cds.all.fasta"
        for ref in sgp grp ogp; do
            cat "${taskdir}/${ref}.cds.all.fasta" >> "${taskdir}/reference.cds.all.fasta"
        done
        hits_sgp=($(cut -f 1 "${taskdir}/hits_slice.tsv" | sort -u))
        python3 -m fasp slice_records_by_ids \
            "${taskdir}/sgp.cds.all.fasta" \
            "${taskdir}/sgp.cds.fasta" \
            "${hits_sgp[@]}"
        hits_reference=($(cut -f 2 "${taskdir}/hits_slice.tsv" | sort -u))
        python3 -m fasp slice_records_by_ids \
            "${taskdir}/reference.cds.all.fasta" \
            "${taskdir}/reference.cds.fasta" \
            "${hits_reference[@]}"
        makeblastdb \
            -in "${taskdir}/reference.cds.fasta" \
            -out "${taskdir}/database_2" \
            -dbtype nucl -hash_index -parse_seqids
    }

    function blastn_hits() {
        blastn -outfmt 6 -evalue 10 \
            -db "${taskdir}/database_2" \
            -query "${taskdir}/sgp.cds.fasta" \
            -out "${taskdir}/hits_2.tsv"
    }

    function slice_hits_2() {
        python3 -m biotp slice_hits_by_crossover_group \
            "${taskdir}/hits_2.tsv" \
            "${taskdir}/hits_slice_2.tsv"
        python3 -m biotp output_besthit_for_subgroup \
            "${taskdir}/hits_2.tsv" \
            "${taskdir}/besthits_2.tsv"
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
        merge_sgp_fasta_2
        merge_grp_fasta_2
        merge_ogp_fasta_2
        makeblastdb_reference
        blastn_hits
        slice_hits_2
    }

    main "$@"
}


