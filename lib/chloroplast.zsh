#!/usr/bin/env zsh
# Last updated: 2024-10-30

#### Function library for Chloroplast genome
: << 'FUNCTIONS'
estimate_mltree_with_whole_chloroplast_genome: Estimate species trees with multiple sequence alignments of whole genome of chloroplast.
FUNCTIONS


function estimate_mltree_with_whole_chloroplast_genome() {
    function usage() {
        cat <<EOS
Usage: construct_chloroplast_dna_msa <arg1> <arg2> <arg3>

    arg1: num       <- Number of organisms.      (e.g., 5)
    arg2: org       <- Names  of organisms.      (e.g., "Cuscuta japonica" "Cuscuta reflexa" ... "Cuscuta campestris")
    arg3: threads   <- Number of threads to use. (e.g., 2)

EOS
    }

    function parse_args() {
        if [[ $# -lt 5 ]]; then
            usage
        fi

        num_orgs=$1; shift
        orgs=()
        for ((i=1; i<=num_orgs; i++)); do
            orgs[i]="${1// /_}"; shift
        done
        threads=$1
    }

    function merge_fasta() {
        touch ${taskdir}/chloroplast.fasta
        for org in "${orgs[@]}"; do
            genus=${org%%_*}
            tmpfile=$(mktemp)
            python3 -m fasp rename_header \
                ${DATA}/${genus}/${org}.dna.chloroplast.fasta \
                $tmpfile \
                $org \
                ""
            cat $tmpfile >> ${taskdir}/chloroplast.fasta
        done
    }

    function construct_msa() {
        mafft --auto ${taskdir}/chloroplast.fasta > ${taskdir}/chloroplast.aln.fasta
        mafft --auto --clustalout ${taskdir}/chloroplast.fasta > ${taskdir}/chloroplast.aln.clustal
    }

    function estimate_tree() {
        raxml-ng \
            --msa ${taskdir}/chloroplast.aln.fasta \
            --all \
            --model GTR+G+I \
            --bs-trees 500 \
            --threads $threads
    }

    function main() {
        parse_args "$@"
        make_taskdir
        merge_fasta
        construct_msa
        estimate_tree
    }

    main "$@"
}