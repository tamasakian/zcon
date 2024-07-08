#!/usr/bin/env zsh

#### Function libs with RAxML-NG

function construct_mltree_msa_pep_blastp_genus_symbol() {
    ## This is depend on construct_msa_pep_blastp_genus_symbol, mafft.sh
    function usage() {
        cat <<EOS
Usage:  construct_mltree_msa_pep_blastp_genus_symbol <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: pid
    arg5: pnm
    ...

EOS
        exit 1
    }

    function construct_mltree() {
        raxml-ng \
            --msa "${taskdir}/${genus}.${symbol}.pep.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 1000 \
            --threads 6 \
            --redo
    }

    function main() {
        construct_msa_pep_blastp_genus_symbol "$@"
        construct_mltree
    }

    main "$@"
}

function construct_mltree_trimmed_msa_pep_blastp_genus_symbol() {
    ## This is depend on trim_msa_pep_blastp_genus_symbol, trimal.sh
    function usage() {
        cat <<EOS
Usage:  construct_mltree_trimmed_msa_pep_blastp_genus_symbol <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: pid
    arg5: pnm
    ...

EOS
        exit 1
    }

    function construct_mltree() {
        raxml-ng \
            --msa "${taskdir}/${genus}.${symbol}.pep.trim.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 1000 \
            --threads 6 \
            --redo
    }

    function main() {
        trim_msa_pep_blastp_genus_symbol "$@"
        construct_mltree
    }

    main "$@"
}

