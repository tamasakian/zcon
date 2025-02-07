#!/usr/bin/env zsh
# Last updated: 2024-12-03
# Tools: BLAST+
# Function libs with BLAST+
: << 'FUNCTIONS'
protein_blast: Protein-Protein BLAST
nucleotide_blast: Nucleotide-Nucleotide BLAST
FUNCTIONS

function protein_blast() {
    function usage() {
        cat <<EOS
Usage: protein_blast <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7>

    arg1: num_qry   <- Number of query organisms.                   (e.g., 1)
    arg2: qry       <- Names  of query organisms.                   (e.g., "Arabidopsis thaliana")
    arg3: num_ref   <- Number of refer organisms.                   (e.g., 2)
    arg4: ref       <- Namea  of refer organisms.                   (e.g., "Cuscuta australis" "Cuscuta campestris")
    arg5: num_pep   <- Number of query proteins.                    (e.g., 2)
    arg6: peps      <- Accession IDs and Names of query proteins.   (e.g., "NP_001326314.1" "AtFLD" "NP_195315.3" "AtFD")
    arg7: evalue    <- Expectation value threshold for saving hits. (e.g., 1e-5)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 7 ]]; then
            usage
        fi

        num_qry=$1; shift
        qrys=()
        echo "Query organisms:${num_qry}"
        for ((i=1; i<=num_qry; i++)); do
            qrys[i]="${1// /_}"; shift
            echo "${i}:${qrys[i]}"
        done

        num_ref=$1; shift
        refs=()
        echo "Refer organisms:${num_ref}"
        for ((i=1; i<=num_ref; i++)); do
            refs[i]="${1// /_}"; shift
            echo "${i}:${refs[i]}"
        done

        num_pep=$1; shift
        typeset -g -A peps
        echo "Query proteins:${num_pep}"
        for ((i=1; i<=num_pep; i++)); do
            key="$1"; shift
            value="$1"; shift
            peps[$key]="$value"
            echo "${i}:${key}"
            echo "${i}:${peps[$key]}"
        done

        ## default evalue
        evalue=10

        if [[ $# != 0 ]]; then
            evalue=$1; shift
        fi

    }

    function generate_qry_fasta() {
        touch "${taskdir}/qry.all.fasta"
        touch "${taskdir}/qry.fasta"
        for qry in "${qrys[@]}"; do
            genus=${qry%%_*}
            cat "${DATA}/${genus}/${qry}.pep.all.fasta" >> "${taskdir}/qry.all.fasta"
        done

        makeblastdb \
            -in "${taskdir}/qry.all.fasta" \
            -dbtype prot \
            -hash_index \
            -parse_seqids

        for pep_key in ${(@k)peps}; do
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry "$pep_key" \
                -db "${taskdir}/qry.all.fasta" \
                -out "$tmpfile"
            python3 -m fasp rename_header \
                "$tmpfile" \
                "$tmpfile" \
                "${peps[$pep_key]}" \
                ""
            cat "$tmpfile" >> "${taskdir}/qry.fasta"
            rm "$tmpfile"
        done
    }

    function generate_ref_fasta() {
        touch "${taskdir}/ref.all.fasta"
        for ref in "${refs[@]}"; do
            genus=${ref%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                "${DATA}/${genus}/${ref}.pep.all.fasta" \
                "$tmpfile" \
                "$ref"
            cat "$tmpfile" >> "${taskdir}/ref.all.fasta"
            rm "$tmpfile"
        done

        makeblastdb \
            -in "${taskdir}/ref.all.fasta" \
            -dbtype prot \
            -hash_index \
            -parse_seqids
    }

    function main() {
        parse_args "$@"
        make_taskdir
        generate_qry_fasta
        generate_ref_fasta
        blastp \
            -outfmt 7 \
            -evalue $evalue \
            -db "${taskdir}/ref.all.fasta" \
            -query "${taskdir}/qry.fasta" \
            -out "${taskdir}/hits.txt"
    }

    main "$@"
}

function nucleotide_blast() {
    function usage() {
        cat <<EOS
Usage: nucleotide_blast <arg1> <arg2> <arg3> <arg4> <arg5> <arg6>

    arg1: num_qry   <- Number of query organisms.                   (e.g., 1)
    arg2: qry       <- Names  of query organisms.                   (e.g., "Arabidopsis thaliana")
    arg3: num_ref   <- Number of refer organisms.                   (e.g., 2)
    arg4: ref       <- Namea  of refer organisms.                   (e.g., "Cuscuta australis" "Cuscuta campestris")
    arg5: num_pep   <- Number of query proteins.                    (e.g., 2)
    arg6: peps      <- Accession IDs and Names of query proteins.   (e.g., "NP_001326314.1" "AtFLD" "NP_195315.3" "AtFD")
    arg7: evalue    <- Expectation value threshold for saving hits. (e.g., 1e-5)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 7 ]]; then
            usage
        fi

        num_qry=$1; shift
        qrys=()
        echo "Query organisms:${num_qry}"
        for ((i=1; i<=num_qry; i++)); do
            qrys[i]="${1// /_}"; shift
            echo "${i}:${qrys[i]}"
        done

        num_ref=$1; shift
        refs=()
        echo "Refer organisms:${num_ref}"
        for ((i=1; i<=num_ref; i++)); do
            refs[i]="${1// /_}"; shift
            echo "${i}:${refs[i]}"
        done

        num_pep=$1; shift
        typeset -g -A peps
        echo "Query peps:${num_pep}"
        for ((i=1; i<=num_pep; i++)); do
            key="$1"; shift
            value="$1"; shift
            peps[$key]="$value"
            echo "${i}:${key}"
            echo "${i}:${peps[$key]}"
        done

        ## default evalue
        evalue=10

        if [[ $# != 0 ]]; then
            evalue=$1; shift
        fi

    }

    function generate_qry_fasta() {
        touch "${taskdir}/qry.all.fasta"
        touch "${taskdir}/qry.fasta"
        for qry in "${qrys[@]}"; do
            genus=${qry%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp rename_headers_feature \
                "${DATA}/${genus}/${qry}.cds.all.fasta" \
                "$tmpfile" \
                "protein_id"
            cat "$tmpfile" >> "${taskdir}/qry.all.fasta"
        done

        makeblastdb \
            -in "${taskdir}/qry.all.fasta" \
            -dbtype nucl \
            -hash_index \
            -parse_seqids

        for pep_key in ${(@k)peps}; do
            tmpfile=$(mktemp)
            blastdbcmd \
                -entry "$pep_key" \
                -db "${taskdir}/qry.all.fasta" \
                -out "$tmpfile"
            python3 -m fasp rename_header \
                "$tmpfile" \
                "$tmpfile" \
                "${peps[$pep_key]}" \
                ""
            cat "$tmpfile" >> "${taskdir}/qry.fasta"
            rm "$tmpfile"
        done
    }

    function generate_ref_fasta() {
        touch "${taskdir}/ref.all.fasta"
        for ref in "${refs[@]}"; do
            genus=${ref%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp rename_headers_feature \
                "${DATA}/${genus}/${ref}.cds.all.fasta" \
                "$tmpfile" \
                "protein_id"
            python3 -m fasp prefix_to_sequence_ids \
                "$tmpfile" \
                "$tmpfile" \
                "$ref"
            cat "$tmpfile" >> "${taskdir}/ref.all.fasta"
            rm "$tmpfile"
        done

        makeblastdb \
            -in "${taskdir}/ref.all.fasta" \
            -dbtype nucl \
            -hash_index \
            -parse_seqids
    }

    function main() {
        parse_args "$@"
        make_taskdir
        generate_qry_fasta
        generate_ref_fasta
        blastn \
            -outfmt 7 \
            -evalue $evalue \
            -db "${taskdir}/ref.all.fasta" \
            -query "${taskdir}/qry.fasta" \
            -out "${taskdir}/hits.txt"
    }

    main "$@"
}
