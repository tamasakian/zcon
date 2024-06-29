#!/bin/zsh

#### Function libs with MAFFT

function construct_msa_pep_blastp_genus_symbol() {
    ## This function is depend on construct_pep_blastp_genus_symbol, blast.sh
    function usage() {
        cat <<EOS
Usage:  construct_msa_pep_blastp_genus_symbol <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: pid
    arg5: pnm
    ...

EOS
        exit 1
    }

    function construct_msa() {
        mafft "${taskdir}/${genus}.${symbol}.pep.fasta" > "${taskdir}/${genus}.${symbol}.pep.aln"
    }

    function main() {
        construct_pep_blastp_genus_symbol "$@"
        construct_msa
    }

    main "$@"
}



