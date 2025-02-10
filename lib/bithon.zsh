#!/usr/bin/env zsh
# Last updated: 2025-2-7
# Tools: bithon BLAST+ diamond
# Function libs with bithon
: << 'FUNCTIONS'
build_longest_isoform: Build the longest isoform from NCBI/ENSEMBL genome data using bithon.
a2a_blastp: All to All Protein-Protein BLAST.
a2a_diamond: All to All Protein-Protein BLAST using diamond.
FUNCTIONS

function build_longest_isoform() {
    function usage() {
        cat << EOS
Usage: build_longest_isoform <arg1> <arg2>

    arg1: sp_num    <- Number of species.   (e.g., 2)
    arg2: sp_name   <- Name of species.     (e.g., "Cuscuta australis" "Cuscuta campestris")

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 2 ]]; then
            usage
        fi
        sp_num=$1
        sp_names=()
        echo "Species:${sp_num}"
        for ((i=1; i<=sp_num; i++)); do
            sp_names[i]="${2// /_}"; shift
            echo "${i}:${sp_names[i]}"
        done
    }

    function build_fasta() {
        mkdir "${taskdir}/input"
        for sp_name in "${sp_names[@]}"; do
            ## ENSEMBL
            if [[ -e "${DATA}/Ensembl/${sp_name}.pep.all.fasta" ]]; then
                bithon ensgls -i "${DATA}/Ensembl/${sp_name}.pep.all.fasta" -o "${taskdir}/input/${sp_name}.fasta" --header transcript
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
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                "${taskdir}/input/${sp_name}.fasta" \
                "$tmpfile" \
                "$sp_name"
            cat "$tmpfile" >> "${taskdir}/input.fasta"
        done
    }

    function main() {
        parse_args "$@"
        make_taskdir
        build_fasta
        setup_fasta
    }

    main "$@"
}

function a2a_blastp() {
    function usage() {
        cat << EOS
Usage: a2a_blastp <arg1> <arg2> <arg3>

    arg1: sp_name_1 <- Name of species. (e.g., "Cuscuta australis")
    arg2: sp_name_2 <- Name of species. (e.g., "Cuscuta campestris")
    arg3: evalue    <- Expectation value threshold for saving hits. (e.g., 1e-5)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 2 ]]; then
            usage
        fi
        sp_names=()
        for ((i=1; i<=2; i++)); do
            sp_names[i]="${1// /_}"; shift
            echo "${i}:${sp_names[i]}"
        done
        ## default evalue
        evalue=10

        if [[ $# != 0 ]]; then
            evalue=$1; shift
        fi
    }

    function build_fasta() {
        mkdir "${taskdir}/input"
        for sp_name in "${sp_names[@]}"; do
            ## ENSEMBL
            if [[ -e "${DATA}/Ensembl/${sp_name}.pep.all.fasta" ]]; then
                bithon ensgls -i "${DATA}/Ensembl/${sp_name}.pep.all.fasta" -o "${taskdir}/input/${sp_name}.fasta" --header transcript
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
            sp_part=(${(s:_:)sp_name})
            sp="${sp_part[1]:0:1}${sp_part[2]}"
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                "${taskdir}/input/${sp_name}.fasta" \
                "$tmpfile" \
                "$sp"
            cat "$tmpfile" >> "${taskdir}/input.fasta"
        done

        makeblastdb \
            -in "${taskdir}/input.fasta" \
            -dbtype prot \
            -hash_index \
            -parse_seqids
    }

    function main() {
        parse_args "$@"
        make_taskdir
        build_fasta
        setup_fasta
        blastp \
            -outfmt 7 \
            -evalue $evalue \
            -db "${taskdir}/input.fasta" \
            -query "${taskdir}/input.fasta" \
            -out "${taskdir}/hits.txt"
    }

    main "$@"
}

function a2a_diamond() {
    function usage() {
        cat << EOS
Usage: a2a_diamond <arg1> <arg2> <arg3>

    arg1: sp_name_1 <- Name of species. (e.g., "Cuscuta australis")
    arg2: sp_name_2 <- Name of species. (e.g., "Cuscuta campestris")
    arg3: max_ts    <- Max target sequnces. (e.g., 10)

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 2 ]]; then
            usage
        fi
        sp_names=()
        for ((i=1; i<=2; i++)); do
            sp_names[i]="${1// /_}"; shift
            echo "${i}:${sp_names[i]}"
        done
        ## default evalue
        max_ts=10

        if [[ $# != 0 ]]; then
            max_ts=$1; shift
        fi
    }

    function build_fasta() {
        mkdir "${taskdir}/input"
        for sp_name in "${sp_names[@]}"; do
            ## ENSEMBL
            if [[ -e "${DATA}/Ensembl/${sp_name}.pep.all.fasta" ]]; then
                bithon ensgls -i "${DATA}/Ensembl/${sp_name}.pep.all.fasta" -o "${taskdir}/input/${sp_name}.fasta" --header transcript
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
            sp_part=(${(s:_:)sp_name})
            sp="${sp_part[1]:0:1}${sp_part[2]}"
            tmpfile=$(mktemp)
            python3 -m fasp prefix_to_sequence_ids \
                "${taskdir}/input/${sp_name}.fasta" \
                "$tmpfile" \
                "$sp"
            cat "$tmpfile" >> "${taskdir}/input.fasta"
        done

        diamond makedb \
            --in ${taskdir}/input.fasta \
            --db ${taskdir}/database
    }

    function main() {
        parse_args "$@"
        make_taskdir
        build_fasta
        setup_fasta
        diamond blastp \
            --db ${taskdir}/database \
            --query ${taskdir}/input.fasta \
            --out ${taskdir}/hits.txt \
            --max-target-seqs $max_ts
    }

    main "$@"
}