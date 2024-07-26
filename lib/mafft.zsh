#!/usr/bin/env zsh

#### Function libs with MAFFT

function construct_msa_pep_blastp_genus_symbol() {
    ## This function is based on construct_pep_blastp_genus_symbol
    ## which is located in blast.zsh
    function usage() {
        cat <<EOS
Usage:  construct_msa_pep_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> (<arg5> <arg6> ...)

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

    function construct_msa() {
        mafft "${taskdir}/${genus}.${symbol}.pep.fasta" > "${taskdir}/${genus}.${symbol}.pep.aln"
    }

    function main() {
        construct_pep_blastp_genus_symbol "$@"
        construct_msa
    }

    main "$@"
}

function addn_construct_msa_pep_blastp_genus_symbol() {
    ## This function is based on addn_construct_pep_blastp_genus_symbol
    ## which is located in blast.zsh
    function usage() {
        cat <<EOS
Usage:  addn_construct_msa_pep_blastp_genus_symbol <arg1> <arg2> (<arg3> <arg4> ... <arg(n+2)>) <arg(n+3)> <arg(n+4)> <arg(n+5)> (<arg(n+6)> <arg(n+7)> ...)

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

    function construct_msa() {
        mafft "${taskdir}/${cons_genus}.${symbol}.pep.fasta" > "${taskdir}/${cons_genus}.${symbol}.pep.aln"
    }

    function main() {
        addn_construct_pep_blastp_genus_symbol "$@"
        construct_msa
    }

    main "$@" 
}

function construct_msa_upstream_region_blastp_genus_symbol() {
    ## This function is based on construct_dna_upstream_region_blastp_genus_symbol
    ## which is located in blast.zsh
    function usage() {
        cat <<EOS
Usage:  construct_msa_upstream_region_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> (<arg6> <arg7> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue
    arg5: bp

    arg6: protein_id
    arg7: protein_name
    ...
    

EOS
        exit 1
    }

    function construct_msa() {
        mafft "${taskdir}/${genus}.${symbol}.dna_upstream_region.fasta" > "${taskdir}/${genus}.${symbol}.dna_upstream_region.aln"
    }

    function main() {
        construct_dna_upstream_region_blastp_genus_symbol "$@"
        construct_msa
    }

    main "$@"
}

function construct_msa_intron_blastp_genus_symbol() {
    ## This function is based on construct_dna_intron_blastp_genus_symbol
    ## which is located in blast.zsh
    function usage() {
        cat <<EOS
Usage:  construct_msa_intron_blastp_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> (<arg6> <arg7> ...)

    arg1: genus
    arg2: symbol
    arg3: symbol_org
    arg4: evalue
    arg5: num

    arg5: protein_id
    arg6: protein_name
    ...
    
EOS
        exit 1
    }

    function construct_msa() {
        for ((i=1; i<=$num; i++)); do
            mafft "${taskdir}/${genus}.intron.${i}.fasta" > "${taskdir}/${genus}.intron.${i}.aln"
        done 
    }

    function main() {
        construct_dna_intron_blastp_genus_symbol "$@"
        construct_msa
    }

    main "$@"
}

