#!/usr/bin/env zsh

#### Function libs with MAFFT

function construct_msa_pep_blastp_genus_symbol() {
    ## This function is based on construct_pep_blastp_genus_symbol, blast.zsh
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

function addn_contruct_msa_pep_blastp_genus_symbol() {
    ## This function is based on addn_construct_pep_blastp_genus_symbol, blast.zsh
    function usage() {
        cat <<EOS
Usage:  addn_contruct_msa_pep_blastp_genus_symbol <arg1> <arg2> (<arg3> <arg4> <arg5> ...) <arg1> <arg2> <arg3> (<arg4> <arg5> ...)

    arg1: cons_genus
    arg2: num [number of add_genus]
    arg3: add1_genus
    arg4: add2_genus
    arg5: add3_genus
    ...
    arg1: symbol
    arg2: symbol_org
    arg3: evalue
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
        addn_construct_pep_blastp_genus_symbol
        construct_msa
    }

    main "$@" 
}


