#!/usr/bin/env zsh

#### Function libs with trimAl

function trim_msa_pep_blastp_genus_symbol() {
    ## This is depend on construct_msa_pep_blastp_genus_symbol
    ## which is located in mafft.zsh
    function usage() {
        cat <<EOS
Usage:  trim_msa_pep_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> (<arg5> <arg6> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue

    arg5: pid
    arg6: pnm
    ...

EOS
        exit 1
    }

    function trim_msa() {
        trimal \
            -in "${taskdir}/${genus}.${symbol}.pep.aln" \
            -out "${taskdir}/${genus}.${symbol}.pep.trim.aln" \
            -automted1
    }

    function main() {
        construct_msa_pep_blastp_genus_symbol "$@"
        trim_msa
    }

    main "$@"
}

function addn_trim_msa_pep_blastp_genus_symbol() {
    ## This function is based on addn_construct_msa_pep_blastp_genus_symbol
    ## which is located in mafft.zsh
    function usage() {
        cat <<EOS
Usage:  addn_trim_msa_pep_blastp_genus_symbol <arg1> <arg2> (<arg3> <arg4> ... <arg(n+2)>) <arg(n+3)> <arg(n+4)> <arg(n+5)> (<arg(n+6)> <arg(n+7)> ...)

    arg1: cons_genus
    arg2: num [number of add_genus]
    arg3: add1_genus
    arg4: add2_genus
    arg(n+2): addn_genus

    arg(n+3): symbol
    arg(n+4): symbol_org
    arg(n+5): evalue

    arg(n+6): pid
    arg(n+7): pnm
    ...

EOS
        exit 1
    }

    function trim_msa() {
        trimal \
            -in "${taskdir}/${cons_genus}.${symbol}.pep.aln" \
            -out "${taskdir}/${cons_genus}.${symbol}.pep.trim.aln" \
            -automted1
    }

    function main() {
        addn_construct_msa_pep_blastp_genus_symbol "$@"
        trim_msa
    }

    main "$@"
}

function trim_sco_msa() {
    ## This function is based on construct_sco_msa
    ## which is located in mafft.zsh
    function usage() {
        cat <<EOS
Usage: trim_sco_msa <arg1> <arg2> ...

    arg1: num (more than 3)
    arg2: org
    ...

EOS
        exit 1
    }

    function trim_msa() {
        scodir="${taskdir}/OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences"
        for file in ${scodir}/*.aln; do
            filename=${file:t:r}
            trimal \
                -in "${scodir}/${filename}.aln" \
                -out "${scodir}/${filename}.trim.aln" \
                -automted1
        done
    }

    function main() {
        construct_sco_msa "$@"
        trim_msa
    }

    main "$@"
}