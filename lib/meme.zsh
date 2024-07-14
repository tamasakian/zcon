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

function fimo_dna_upstream_region_genus_symbol_test() {
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

