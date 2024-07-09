#!/usr/bin/env zsh

#### Function libs with datasets

function fetch_genome_by_genus() {
    function usage() {
        cat <<EOS
Usage:  fetch_genome_by_genus <arg1>

    arg1: genus

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 1 ]]; then
            usage
        fi
        genus="$1"
    }

    function make_dir() {
        mkdir -p "${DOWNLOAD}/${genus}"
        mkdir -p "${DATA}/${genus}"
    } 

    function download_genome_by_genus() {
        datasets download genome taxon "$genus" \
            --include genome,protein,cds,gff3 \
            --reference --annotated \
            --filename "${DOWNLOAD}/${genus}/dataset.zip"
        unzip -o "${DOWNLOAD}/${genus}/dataset.zip" -d "${DOWNLOAD}/${genus}"
    }
    
    function declare_genome_by_genus() {
        acc_li=() org_li=() asm_li=()

        tmpfile=$(mktemp)
        python3 -m biotp output_acc_org_asm \
            "${DOWNLOAD}/${genus}/ncbi_dataset/data/assembly_data_report.jsonl" \
            > "$tmpfile"

        while IFS=" " read -r i acc org asm; do
            acc_li[i]="$acc" # acc: WGS accession
            org_li[i]="$org" # org: Scientific name
            asm_li[i]="$asm" # asm: Assembly
            echo "${i}, acc: ${acc_li[i]}, org: ${org_li[i]}, asm: ${asm_li[i]}"
        done < "$tmpfile"

        rm "$tmpfile"
    }

    function send_genome_to_data() {
        for i in "${(k)org_li}"; do
            echo "$i"
            acc=${acc_li[$i]}; org=${org_li[$i]}; asm=${asm_li[$i]}
            cp "${DOWNLOAD}/${genus}/ncbi_dataset/data/${acc}/${acc}_${asm}_genomic.fna" \
                "${DATA}/${genus}/${org}.dna.toplevel.fasta"
            cp "${DOWNLOAD}/${genus}/ncbi_dataset/data/${acc}/cds_from_genomic.fna" \
                "${DATA}/${genus}/${org}.cds.all.fasta"
            cp "${DOWNLOAD}/${genus}/ncbi_dataset/data/${acc}/protein.faa" \
                "${DATA}/${genus}/${org}.pep.all.fasta"
            cp "${DOWNLOAD}/${genus}/ncbi_dataset/data/${acc}/genomic.gff" \
                "${DATA}/${genus}/${org}.genome.gff"
        done
    }

    function main() {
        parse_args "$@"
        make_dir
        download_genome_by_genus
        declare_genome_by_genus
        send_genome_to_data
    }

    main "$@"
}

function redeclare_genome_by_genus() {
    function usage() {
        cat <<EOS
Usage: redeclare_genome_by_genus <arg1>

    arg1: genus

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 1 ]]; then
            usage
        fi
        genus="$1"
    }

    function declare_genome_by_genus() {
        acc_li=() org_li=() asm_li=()
        tmpfile=$(mktemp)
        python3 -m biotp output_acc_org_asm \
            "${DOWNLOAD}/${genus}/ncbi_dataset/data/assembly_data_report.jsonl" \
            > "$tmpfile"
        while IFS=" " read -r i acc org asm; do
            acc_li[i]="$acc" 
            org_li[i]="$org" 
            asm_li[i]="$asm"
            echo "${i}, acc: ${acc_li[i]}, org: ${org_li[i]}, asm: ${asm_li[i]}"
        done < "$tmpfile"
        rm "$tmpfile"
    }

    function main() {
        parse_args "$@"
        declare_genome_by_genus
    }

    main "$@"
}

function fetch_gene_by_symbol() {
    function usage() {
        cat <<EOS
Usage:  fetch_gene_by_symbol <arg1> <arg2>

    arg1: symbol
    arg2: org

EOS
        exit 1
    }

    function parse_args() {
        if [[ $# != 2 ]]; then
            usage
        fi
        symbol="$1"
        org="$2"
        org_us="${2/ /_}"
    }

    function make_dir() {
        mkdir -p "${DOWNLOAD}/${org_us}"
        mkdir -p "${DOWNLOAD}/${org_us}/${symbol}"
        mkdir -p "${DATA}/${org_us}"
    }

    function download_gene_by_symbol() {
        datasets download gene symbol "$symbol" \
            --taxon "$org" \
            --include protein,cds \
            --filename "${DOWNLOAD}/${org_us}/${symbol}/dataset.zip"
        unzip -o "${DOWNLOAD}/${org_us}/${symbol}/dataset.zip" -d "${DOWNLOAD}/${org_us}/${symbol}"
    }

    function send_gene_to_data() {
        cp "${DOWNLOAD}/${org_us}/${symbol}/ncbi_dataset/data/cds.fna" \
            "${DATA}/${org_us}/${symbol}.cds.fasta"
        cp "${DOWNLOAD}/${org_us}/${symbol}/ncbi_dataset/data/protein.faa" \
            "${DATA}/${org_us}/${symbol}.pep.fasta"
    }

    function main() {
        parse_args "$@"
        make_dir
        download_gene_by_symbol
        send_gene_to_data
    }

    main "$@"
}