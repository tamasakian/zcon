#!/bin/zsh

#### Function libs with trimAl

function trim_msa_pep_blastp_genus_symbol() {
    ## This is depend on construct_msa_pep_blastp_genus_symbol, mafft.sh
    function usage() {
        cat <<EOS
Usage:  trim_msa_pep_blastp_genus_symbol <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: pid
    arg5: pnm
    ...

EOS
        exit 1
    }

    function trim_msa() {
        trimal \
            -in "${taskdir}/${genus}.${symbol}.pep.aln" \
            -out "${taskdir}/${genus}.${symbol}.pep.trim.aln"
            -automted1
    }

    function main() {
        construct_msa_pep_blastp_genus_symbol "$@"
        trim_msa
    }

    main "$@"
}