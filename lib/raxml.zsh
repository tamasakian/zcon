#!/usr/bin/env zsh

#### Function libs with RAxML-NG

function construct_mltree_msa_pep_blastp_genus_symbol() {
    ## This is depend on construct_msa_pep_blastp_genus_symbol
    ## which is located in mafft.zsh
    function usage() {
        cat <<EOS
Usage:  construct_mltree_msa_pep_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> (<arg5> <arg6> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue

    arg5: protein_id
    arg6: protein_name
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
    ## This is depend on trim_msa_pep_blastp_genus_symbol
    ## which is located in trimal.zsh
    function usage() {
        cat <<EOS
Usage:  construct_mltree_trimmed_msa_pep_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> (<arg5> <arg6> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue

    arg5: protein_id
    arg6: protein_name
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

function addn_construct_mltree_msa_pep_blastp_genus_symbol() {
    ## This function is based on addn_construct_msa_pep_blastp_genus_symbol
    ## which is located in mafft.zsh
    function usage() {
        cat <<EOS
Usage:  addn_construct_mltree_msa_pep_blastp_genus_symbol <arg1> <arg2> (<arg3> <arg4> ... <arg(n+2)>) <arg(n+3)> <arg(n+4)> <arg(n+5)> (<arg(n+6)> <arg(n+7)> ...)

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

    function construct_mltree() {
        raxml-ng \
            --msa "${taskdir}/${cons_genus}.${symbol}.pep.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 1000 \
            --threads 6 \
            --redo
    }

    function main() {
        addn_construct_msa_pep_blastp_genus_symbol "$@"
        construct_mltree
    }

    main "$@"
}

function addn_construct_mltree_trimmed_msa_pep_blastp_genus_symbol() {
    ## This function is based on addn_trim_msa_pep_blastp_genus_symbol
    ## which is located in trimal.zsh
    function usage() {
        cat <<EOS
Usage:  addn_construct_mltree_trimmed_msa_pep_blastp_genus_symbol <arg1> <arg2> (<arg3> <arg4> ... <arg(n+2)>) <arg(n+3)> <arg(n+4)> <arg(n+5)> (<arg(n+6)> <arg(n+7)> ...)

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

    function construct_mltree() {
        raxml-ng \
            --msa "${taskdir}/${cons_genus}.${symbol}.pep.trim.aln" \
            --all \
            --model Blosum62 \
            --bs-trees 1000 \
            --threads 6 \
            --redo
    }

    function main() {
        addn_trim_msa_pep_blastp_genus_symbol "$@"
        construct_mltree
    }

    main "$@"
}
