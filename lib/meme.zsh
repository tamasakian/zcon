#!/usr/bin/env zsh

#### Function libs with MEME

function download_meme_databases() {

    function download_databases() {
        wget -P "$DOWNLOAD" "https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.24.tgz"
        wget -P "$DOWNLOAD" "https://meme-suite.org/meme/meme-software/Databases/gomo/gomo_databases.3.2.tgz"
        wget -P "$DOWNLOAD" "https://meme-suite.org/meme/meme-software/Databases/tgene/tgene_databases.1.0.tgz"
    }

    function copy_databases() {
        tar -zxvf "${DOWNLOAD}/motif_databases.12.24.tgz" -C "$DATA"
        tar -zxvf "${DOWNLOAD}/gomo_databases.3.2.tgz"    -C "$DATA"
        tar -zxvf "${DOWNLOAD}/tgene_databases.1.0.tgz"   -C "$DATA"
    }

    function main() {
        download_databases
        copy_databases
    }
    
    main
}

function fimo_dna_upstream_region() {
    ## This function is based on construct_dna_upstream_region_genus_symbol
    ## which is located in blast.zsh
    function usage() {
        cat <<EOS
Usage:  fimo_dna_upstream_region <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> <arg7> (<arg8> <arg9> ...)

    arg1: motif_filename
    arg2: qvalue (fimo)

    arg3: genus
    arg4: symbol
    arg5: symbol_org
    arg6: evalue (blastn)
    arg7: bp
    arg8: pid
    arg9: pnm
    ...

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 9 ]]; then
            usage
        fi
        motif_filename="$1"
        qvalue=$2
    }

    function find_individual_motif_occurrences() {
        fimo \
            --o "${taskdir}/fimo_out" \
            --thresh $qvalue \
            --qv-thresh \
            --verbosity 1 \
            "$motif_filename" \
            "${taskdir}/${genus}.${symbol}.dna_upstream_region.fasta"
    }
    
    function main() {
        parse_args "$@"
        construct_dna_upstream_region_genus_symbol "${@:3}"
        find_individual_motif_occurrences
    }

    main "$@"
}

function fimo_protein_genus_symbol() {
    ## This function is based on construct_pep_blastp_genus_symbol
    ## which is located in blast.zsh
    function usage() {
        cat <<EOS
Usage:  fimo_protein_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> (<arg7> <arg8> ...)

    arg1: motif_filename
    arg2: qvalue (fimo)

    arg3: genus
    arg4: symbol
    arg5: symbol_org
    arg6: evalue (blastn)
    arg7: pid
    arg8: pnm
    ...

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 8 ]]; then
            usage
        fi
        motif_filename="$1"
        qvalue=$2
    }

    function find_individual_motif_occurrences() {
        fimo \
            --o "${taskdir}/fimo_out" \
            --thresh $qvalue \
            --qv-thresh \
            --verbosity 1 \
            "$motif_filename" \
            "${taskdir}/${genus}.${symbol}.pep.fasta"
    }
    
    function main() {
        parse_args "$@"
        construct_pep_blastp_genus_symbol "${@:3}"
        find_individual_motif_occurrences
    }

    main "$@"
}

function fimo_addn_protein_genus_symbol() {
    ## This function is based on addn_construct_pep_blastp_genus_symbol
    ## which is located in blast.zsh
    function usage() {
        cat <<EOS
Usage:  fimo_protein_genus_symbol <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> ...

    arg1: motif_filename
    arg2: qvalue (fimo)

    arg3: cons_genus
    arg4: num [number of add_genus]
    arg5: add1_genus
    arg6: add2_genus
    arg(n+4): addn_genus

    arg(n+5): symbol
    arg(n+6): symbol_org
    arg(n+7): evalue

    arg(n+8): protein_id
    arg(n+9): protein_name
    ...

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# -lt 8 ]]; then
            usage
        fi
        motif_filename="$1"
        qvalue=$2
    }

    function find_individual_motif_occurrences() {
        fimo \
            --o "${taskdir}/fimo_out" \
            --thresh $qvalue \
            --qv-thresh \
            --verbosity 1 \
            "$motif_filename" \
            "${taskdir}/${cons_genus}.${symbol}.pep.fasta"
    }
    
    function main() {
        parse_args "$@"
        addn_construct_pep_blastp_genus_symbol "${@:3}"
        find_individual_motif_occurrences
    }

    main "$@"
}
