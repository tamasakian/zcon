#!/usr/bin/env zsh

#### Function library for HGT research
# Last updated: 2024-10-30

function detect_crossover_besthits() {
    function usage() {
        cat <<EOS
Usage: detect_crossover_besthits <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> 

    arg1: num_sgp   <- Number of organisms in subgroup (e.g., 3)
    arg2: org_sgp   <- Names  of organisms in subgroup (e.g., "Cuscuta australis" "Cuscuta campestris" ... "Cuscuta europaea")
    arg3: num_grp   <- Number of organisms in group
    arg4: org_grp   <- Names  of organisms in group    (e.g., "species1" "species2" ... "speciesN")
    arg5: num_ogp   <- Number of organisms in outgroup
    arg6: org_ogp   <- Names  of organisms in outgroup (e.g., "species1" "species2" ... "speciesN")

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 6 ]]; then
            usage
        fi

        num_sgp=$1; shift
        org_sgp=()
        for ((i=1; i<=num_sgp; i++)); do
            org_sgp[i]="${1// /_}"; shift
            echo "${i}: ${org_sgp[i]}"
        done

        num_grp=$1; shift
        org_grp=()
        for ((j=1; j<=num_grp; j++)); do
            org_grp[j]="${1// /_}"; shift
            echo "${j}: ${org_grp[j]}"
        done

        num_ogp=$1; shift
        org_ogp=()
        for ((k=1; k<=num_ogp; k++)); do
            org_ogp[k]="${1// /_}"; shift
            echo "${k}: ${org_ogp[k]}"
        done

        max_target_seqs=$(((num_sgp + num_grp + num_ogp) * 10))
        echo "max_target_seqs: ${max_target_seqs}"
    }

    function merge_sgp_fasta() {
        touch ${taskdir}/sgp.fasta
        for org in "${org_sgp[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                ${DATA}/${genus}/${org}.pep.all.fasta \
                $tmpfile \
                sgp_${org}_
            cat $tmpfile >> ${taskdir}/sgp.fasta
            rm $tmpfile
        done
    }

    function merge_grp_fasta() {
        touch ${taskdir}/grp.fasta
        for org in "${org_grp[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                ${DATA}/${genus}/${org}.pep.all.fasta \
                $tmpfile \
                grp_${org}_
            cat $tmpfile >> ${taskdir}/grp.fasta
            rm $tmpfile
        done
    }

    function merge_ogp_fasta() {
        touch ${taskdir}/ogp.fasta
        for org in "${org_ogp[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                ${DATA}/${genus}/${org}.pep.all.fasta \
                $tmpfile \
                ogp_${org}_
            cat $tmpfile >> ${taskdir}/ogp.fasta
            rm $tmpfile
        done
    }

    function make_database() {
        touch ${taskdir}/reference.fasta
        for ref in sgp grp ogp; do
            cat ${taskdir}/${ref}.fasta >> ${taskdir}/reference.fasta
        done
        diamond makedb \
            --in ${taskdir}/reference.fasta \
            --db ${taskdir}/database
    }

    function search_homology() {
        diamond blastp \
            --db ${taskdir}/database \
            --query ${taskdir}/sgp.fasta \
            --out ${taskdir}/hits.tsv \
            --max-target-seqs $max_target_seqs
    }

    function detect_crossover() {
        python3 -m biotp slice_hits_by_crossover_group \
            ${taskdir}/hits.tsv \
            ${taskdir}/hits_slice.tsv
        python3 -m biotp output_besthit_for_subgroup \
            ${taskdir}/hits.tsv \
            ${taskdir}/besthits.tsv
    }

    function merge_sgp_fasta_cds() {
        touch ${taskdir}/sgp.cds.all.fasta
        for org in "${org_sgp[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m fasp rename_headers_feature \
                ${DATA}/${genus}/${org}.cds.all.fasta \
                $tmpfile1 \
                protein_id
            python3 -m fasp prefix_to_sequence_ids \
                $tmpfile1 \
                $tmpfile2 \
                sgp_${org}_
            cat $tmpfile2 >> ${taskdir}/sgp.cds.all.fasta
            rm $tmpfile1 $tmpfile2
        done
    }

    function merge_grp_fasta_cds() {
        touch ${taskdir}/grp.cds.all.fasta
        for org in "${org_grp[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m fasp rename_headers_feature \
                ${DATA}/${genus}/${org}.cds.all.fasta \
                $tmpfile1 \
                protein_id
            python3 -m fasp prefix_to_sequence_ids \
                $tmpfile1 \
                $tmpfile2 \
                grp_${org}_
            cat $tmpfile2 >> ${taskdir}/grp.cds.all.fasta
            rm $tmpfile1 $tmpfile2
        done
    }

    function merge_ogp_fasta_cds() {
        touch ${taskdir}/ogp.cds.all.fasta
        for org in "${org_ogp[@]}"; do
            genus=${org%%_*}
            tmpfile1=$(mktemp)
            tmpfile2=$(mktemp)
            python3 -m fasp rename_headers_feature \
                ${DATA}/${genus}/${org}.cds.all.fasta \
                $tmpfile1 \
                protein_id
            python3 -m fasp prefix_to_sequence_ids \
                $tmpfile1 \
                $tmpfile2 \
                ogp_${org}_
            cat $tmpfile2 >> ${taskdir}/ogp.cds.all.fasta
            rm $tmpfile1 $tmpfile2
        done
    }

    function make_database_cds() {
        touch ${taskdir}/reference.cds.all.fasta
        for ref in sgp grp ogp; do
            cat ${taskdir}/${ref}.cds.all.fasta >> ${taskdir}/reference.cds.all.fasta
        done
        tmpfile1=$(mktemp)
        tmpfile2=$(mktemp)
        cut -f 1 ${taskdir}/hits_slice.tsv | sort -u | xargs -n 1000 > $tmpfile1
        cut -f 2 ${taskdir}/hits_slice.tsv | sort -u | xargs -n 1000 > $tmpfile2
        python3 -m fasp slice_records_by_idfile \
            ${taskdir}/sgp.cds.all.fasta \
            ${taskdir}/sgp.cds.fasta \
            $tmpfile1
        python3 -m fasp slice_records_by_idfile \
            ${taskdir}/reference.cds.all.fasta \
            ${taskdir}/reference.cds.fasta \
            $tmpfile2
        rm $tmpfile1 $tmpfile2
        python3 -m fasp assign_unique_ids \
            ${taskdir}/reference.cds.fasta \
            ${taskdir}/reference.cds.unique.fasta
        makeblastdb \
            -in ${taskdir}/reference.cds.unique.fasta \
            -out ${taskdir}/database_cds \
            -dbtype nucl -hash_index -parse_seqids
    }

    function search_homology_cds() {
        blastn -outfmt 6 -evalue 10 \
            -db ${taskdir}/database_cds \
            -query ${taskdir}/sgp.cds.fasta \
            -out ${taskdir}/hits_cds.tsv \
            -max_target_seqs $max_target_seqs
    }

    function detect_crossover_cds() {
        python3 -m biotp slice_hits_by_crossover_group \
            ${taskdir}/hits_cds.tsv \
            ${taskdir}/hits_slice_cds.tsv
        python3 -m biotp output_besthit_for_subgroup \
            ${taskdir}/hits_cds.tsv \
            ${taskdir}/besthits_cds.tsv
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_sgp_fasta
        merge_grp_fasta
        merge_ogp_fasta
        make_database
        search_homology
        detect_crossover
        merge_sgp_fasta_cds
        merge_grp_fasta_cds
        merge_ogp_fasta_cds
        make_database_cds
        search_homology_cds
        detect_crossover_cds
    }

    main "$@"
}


